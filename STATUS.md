# Project status

Status date: **2026-08-12**
Current milestone: **M6 arithmetic/oracle baseline with x4 performance decision frozen**

## Executive assessment

The project now has a verification-oriented arithmetic baseline rather than a
blanket replacement for `MultiFloats.jl`:

- x2/x3 fixed-cost direct-FMA networks with pinned structural verification;
- x4 cancellation-safe direct FMA;
- exact/adaptive higher-precision reference oracles;
- correctness-first x5-x8 add/mul/direct-FMA baselines;
- an accepted Float64x5 reciprocal seed/correction path;
- permanent BigFloat128/256/512 versus Float64x2/x4/x8 comparison benchmarks.

The x5-x8 safe mul/FMA implementations are intentionally over-complete reference
implementations. They are not performance kernels.

## x4 performance decision

A fresh Zen 3 diagnostic tested whether the conservative x4 repair could be
recovered cheaply without hidden FastTwoSum assumptions.

### What was tested

1. profile the final four-term `renormalize` input;
2. unroll several renormalization passes with the authoritative fallback;
3. test a scalar guarded tail where each FastTwoSum is used only when its
   `|a| >= |b|` precondition is explicitly true, otherwise TwoSum is used;
4. keep the current safe renormalization fallback whenever the four-term result
   is not already stable;
5. require bitwise equality with the current safe direct-FMA result on ordinary,
   destructive-cancellation, and wide-exponent corpora.

### Result

For the final-renormalization input:

- ordinary: 20,000 / 20,000 cases already stable before generic renormalization;
- destructive cancellation: 2,364 / 5,000 already stable, 2,636 / 5,000 need
  only one changing pass;
- a simple unrolled-renormalization specialization improves scalar time by only
  about 4% and remains about 1.9x the time of upstream `x*y+c`.

The stronger guarded-FastTwoSum candidate passes its empirical safety gates:

- 20,000 ordinary cases: bitwise equal to current safe;
- 10,000 destructive-cancellation cases: zero normalization failures and zero
  guarded-vs-safe mismatches;
- 20,000 wide-exponent cases: zero guarded-vs-safe mismatches.

But it provides **no performance win** on the hosted Zen 3 run:

| Float64x4 scalar, 20k cases | Time |
|---|---:|
| historical unsafe direct tail | 0.156 ms |
| current safe direct FMA | 0.658 ms |
| guarded direct FMA | 0.660 ms |
| upstream `x*y+c` | 0.329 ms |

Therefore the x4 scalar slowdown is not primarily the generic renormalization
loop, and dynamic FastTwoSum guards do not recover it. **Further x4 scalar
micro-optimization is stopped.** A future x4 change is justified only by a new
source-specific fixed-cost proof/compressor, not by more guard or renormalization
tuning.

The SIMD conclusion is different. On the same class of Zen 3 hosted runs,
current direct x4 FMA remains roughly 1.3x-1.7x faster than upstream mul+add for
Vec2/Vec4/Vec8 workloads. This is the useful x4 deployment path.

## x4 numerical comparison with original MultiFloats composition

Upstream `x*y+c` is not being characterized as unstable: it is the normal
separately-rounded MultiFloats multiplication and addition path. The project
provides a different direct-FMA operation.

On 5,000 seeded destructive-cancellation cases, using exact `Rational{BigInt}`
`x*y+c` as the error reference:

- project direct FMA had smaller exact error in **3,215 / 5,000** cases;
- upstream `x*y+c` had smaller error in **591 / 5,000** cases;
- they had equal error in **1,194 / 5,000** cases.

Thus direct FMA is often more accurate under severe cancellation, but it is not
universally more accurate than the separately-rounded composition.

## Current comparison with released MultiFloats v3.2.6

The project currently pins MultiFloats v3.2.6. Ordinary x2/x4 add and multiply
remain upstream MultiFloats strengths; this project should not replace them just
for the sake of replacement.

Recent hosted scalar throughput-equivalent measurements:

| Operation | Original/upstream path | Project path | Decision |
|---|---:|---:|---|
| x2 mul+add | ~2.1 ns/op | `fma_fast` ~1.5 ns/op | project direct FMA faster |
| x4 add | ~7.1 ns/op | upstream is the implementation | keep upstream |
| x4 mul | ~7.3 ns/op | upstream is the implementation | keep upstream |
| x4 mul+add / direct FMA | ~16.3 ns/op | ~32.9 ns/op | upstream faster scalar |
| x4 Vec2/4/8 mul+add / direct FMA | baseline | project ~1.3x-1.7x faster | project useful for SIMD |

The project is therefore **not** currently a blanket faster replacement for
MultiFloats. It is a direct-FMA / verification / higher-width extension with
selected performance wins.

## BigFloat comparison

The permanent benchmark pairs BigFloat128/256/512 with Float64x2/x4/x8. These
pairs intentionally give BigFloat about 20.8% more nominal significand bits.

Latest Zen 3 snapshot:

- BigFloat128 add/mul: ~33.7 / 47.7 ns; original Float64x2: ~1.4 / 1.0 ns;
- BigFloat256 add/mul: ~41.6 / 69.6 ns; original Float64x4: ~7.1 / 7.3 ns;
- BigFloat512 add/mul: ~51.5 / 107 ns; current correctness-only Float64x8:
  ~285 ns / **367,305 ns**.

The x8 number is a property of the current over-complete reference algorithm,
not a reason to conclude that optimized fixed-limb x8 cannot work.

## x8 optimization decision

**Optimize x8, but redesign it; do not micro-optimize the present safe
implementation.**

The current x8 multiplication constructs all 128 TwoProd product/residual
components, orders them, fully renormalizes them, and truncates only afterward.
No local finalizer optimization can close a three-orders-of-magnitude-plus gap.

The performance track is therefore:

1. keep the current x5-x8 exact/safe implementations as independent rejection
   oracles;
2. derive a diagonal/truncated fixed-cost multiplication network, beginning at
   x5/x6 where proof/debug cost is lower;
3. build the matching direct-FMA network without rounded mul-then-add
   composition;
4. target SIMD/`MultiFloatVec` from the start, because x4 evidence shows the
   strongest fixed-limb performance advantage there;
5. promote the design to x8 only after x5/x6 achieve zero oracle failures and a
   material performance win;
6. compare x8 against BigFloat512 and against any relevant upstream master
   implementation before M7 integration.

Go/no-go for production x8: it must cease to be a reference-only kernel and show
an end-to-end material advantage in the intended vector/matrix workload. If it
cannot beat or materially outperform BigFloat512 after the network redesign,
M7 should use x4 or BigFloat rather than forcing x8.

## Current arithmetic decisions

| Component | Decision |
|---|---|
| x2 direct FMA | keep; scalar performance win observed |
| x3 direct FMA | keep; proof-driven repair, mainly SIMD performance case |
| x4 scalar direct FMA | correctness feature; stop micro-optimization |
| x4 SIMD direct FMA | keep; measured performance win |
| x5-x8 `reference_*` | authoritative Experimental rejection oracle |
| x5-x8 safe add/mul/FMA | correctness/reference baselines |
| current x8 safe mul/FMA | never use as M7 hot-path performance kernel |
| optimized x8 | pursue only as a new fixed-cost SIMD-first network |
| rounded mul-safe -> add-safe as direct FMA | rejected under cancellation |
| native DOT/SYRK/TRSM/GEMM | after competitive arithmetic survivors exist |

## Immediate next performance milestone

Do **not** spend another iteration on x4 scalar guards or generic
renormalization. Start the high-limb network-reduction track at x5/x6, preserving
all existing exact/oracle tests. The target is a fixed-cost, allocation-free,
SIMD-friendly multiplication/direct-FMA construction that can later be lifted
to x8 and measured against BigFloat512.
