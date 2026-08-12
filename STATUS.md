# Project status

Status date: **2026-08-12**
Current milestone: **M5 complete baseline — safe direct Float64x5-x8 FMA family**

## Executive assessment

The project now has correctness-first higher-limb baselines for all three core
arithmetic operations needed by later residual/refinement kernels:

- x2/x3 fixed-cost FMA networks with pinned structural verification;
- x4 cancellation-safe TwoSum + `renormalize` FMA baseline;
- exact-rational Experimental add/sub/mul/FMA oracles for x5-x8;
- one no-FastTwoSum safe addition family for Float64x5-x8;
- one no-FastTwoSum safe multiplication family for Float64x5-x8;
- one **direct** safe FMA family for Float64x5-x8 that forms the full exact
  `x*y+c` expansion before any N-limb truncation.

M5 is therefore complete as a correctness baseline. The current x5-x8 safe
mul/FMA implementations are deliberately over-complete, allocating, and slow;
they are reference implementations for later minimization rather than candidate
production kernels.

The requested BigFloat comparison makes that distinction concrete: original
upstream MultiFloat x2/x4 arithmetic is much faster than BigFloat128/256 in the
hosted throughput benchmark, whereas the current correctness-only x8 safe
multiplication is roughly 3,600x slower than BigFloat512. A smaller verified x8
product/FMA network is mandatory before M7 performance claims.

## M5 safe direct FMA family

`Experimental.fma_safe` covers N=5:8. Width N:

1. obtains all `2N^2` exact TwoProd rounded-product/residual components from the
   M4 multiplication decomposition;
2. appends all N exact limbs of `c`;
3. canonicalizes signed zero and deterministically sorts the full component set;
4. fully `renormalize`s the exact direct `x*y+c` expansion;
5. retains the leading N limbs and renormalizes that head.

The full component counts are 55, 78, 105, and 136 for x5, x6, x7, and x8.
No FastTwoSum is used. Product-side overflow/subnormal/underflow exclusions
inherit the current M4 `mul_safe` contract.

The frozen x5 wrapper `fma5_safe` remains available; the generic family is
required to match it bit-for-bit. Width-specific wrappers `fma6_safe`,
`fma7_safe`, and `fma8_safe` call the generic construction.

### Permanent correctness gates

For every accepted width/case CI requires:

```text
value(result) + value(discarded normalized tail)
= value(x) * value(y) + value(c)
```

It also checks every individual TwoProd pair by exact `Rational{BigInt}`
reconstruction, normalized full/output expansions, bitwise `reference_fma`
equality, bitwise x/y symmetry, generic/wrapper equality, identities,
dense/scaled inputs, powers of two, and destructive cancellation.

All x5-x8 direct-FMA gates pass on Linux Julia 1.10/current and macOS current.

### First width-specific direct-FMA diagnostics

Zen 3, Julia 1.10.11. `C=1` remains a deliberately loose empirical regression
gate, not a theorem.

| Width | Max observed `|err|/(u^N(|xy|+|c|))` | Direct oracle / norm / symmetry failures |
|---|---:|---:|
| x5 | 0.0492768 | 0 / 0 / 0 |
| x6 | 0.0227603 | 0 / 0 / 0 |
| x7 | 0.00452122 | 0 / 0 / 0 |
| x8 | 0.000825581 | 0 / 0 / 0 |

The maximum is taken over the ordinary/scaled seeded corpora for each width.
Destructive-cancellation operand-relative constants were many orders smaller
because the scale is `|xy|+|c|` while the surviving result is tiny.

### Direct FMA versus rounded multiplication followed by addition

The exact direct-FMA oracle exposes intermediate-rounding loss even though
`mul_safe` and `add_safe` separately match their own operation-specific oracles:

| Width | Ordinary composition mismatches | Scaled mismatches | Destructive-cancellation mismatches |
|---|---:|---:|---:|
| x5 | 45 / 200 | 41 / 200 | **150 / 150** |
| x6 | 8 / 50 | 14 / 50 | **48 / 48** |
| x7 | 7 / 35 | 9 / 35 | **32 / 32** |
| x8 | 6 / 25 | 6 / 25 | **24 / 24** |

Thus `add_safe(mul_safe(x,y),c)` is permanently rejected as an exact substitute
for direct FMA in cancellation-sensitive residual, refinement, or certificate
work.

### Safe-baseline timing A/B

Same hosted Zen 3 run; case counts differ because the exact expansion grows to
136 components at x8.

| Width | Direct safe FMA | Safe mul-then-add | Composition/direct |
|---|---:|---:|---:|
| x5 | 9.852 ms / 120 | 7.993 ms / 120 | 0.811x |
| x6 | 3.044 ms / 20 | 2.550 ms / 20 | 0.837x |
| x7 | 3.186 ms / 12 | 2.726 ms / 12 | 0.856x |
| x8 | 2.625 ms / 6 | 2.314 ms / 6 | 0.882x |

The direct baseline remains slower, but its relative overhead shrinks with width
in this sample. These timings are reference-path diagnostics, not performance
targets.

## Requested BigFloat comparison

The permanent comparison benchmark pairs BigFloat128/256/512 with Float64x2/x4/
x8. The pairings are intentionally unequal precision: BigFloat has about 20.8%
more nominal significand bits.

First hosted throughput-equivalent results:

- original x2 add/mul were about 36x/91x higher throughput than BigFloat128;
- original x4 add/mul were about 6x/9x higher throughput than BigFloat256;
- current safe x8 add was about 4.6x slower than BigFloat512;
- current safe x8 multiplication was about **3,644x slower** than BigFloat512.

See `docs/BIGFLOAT_MULTIFLOAT_COMPARISON.md`. The x8 result is a development
signal about the over-complete safe implementation, not evidence against an
optimized fixed-limb design.

## M4 safe multiplication family

`Experimental.mul_safe` and `mul5_safe` ... `mul8_safe` retain all `2N^2`
TwoProd components, canonicalize their order, fully renormalize, and truncate
only after the exact product expansion exists. First max observed
`|err|/(u^N|xy|)` values were 0.0484069, 0.0183703, 0.00725543, and 0.00253803
for x5-x8, with zero measured oracle/normalization/commutativity failures.

## M3 safe addition family

`Experimental.add_safe` and `add5_safe` ... `add8_safe` use N general TwoSums,
full exact 2N-term renormalization, N-limb truncation, and head renormalization.
First max observed `|err|/(u^N(|x|+|y|))` values were 0.0528336, 0.0254805,
0.0109878, and 0.00541727 for x5-x8.

## Current decisions

| Component | Decision |
|---|---|
| x2/x3 fixed-cost FMA | structurally verified; architecture-sensitive performance |
| x4 FMA | cancellation-safe correctness baseline |
| x5-x8 `reference_*` | exact Experimental oracle |
| x5-x8 `add_safe` | accepted M3 Experimental correctness family |
| x5-x8 `mul_safe` | accepted M4 Experimental correctness family |
| x5-x8 `fma_safe` | accepted M5 Experimental direct-FMA correctness family |
| mul-safe -> add-safe FMA composition | rejected as direct-FMA substitute; fails every tested destructive-cancellation case |
| safe x8 mul/FMA | correctness reference only; far too slow for M7 hot paths |
| fixed-cost high-limb add/mul/FMA | formal tail analysis + source-specific network search next |
| reciprocal/div/sqrt | M6 correctness baseline next |
| native DOT/SYRK/TRSM/GEMM | M7 after competitive arithmetic survivors exist |

## Immediate next milestone

Start **M6** with a correctness-first Float64x5 reciprocal/division baseline and
an independent high-precision oracle. Use direct FMA/submul for residual
corrections where it avoids intermediate rounding. In parallel, begin the formal
error/tail analysis and product-network reduction needed to replace the current
over-complete x8 multiplication/FMA before M7 performance integration.
