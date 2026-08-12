# Project status

Status date: **2026-08-12**
Current milestone: **M2 — exact x5-x8 reference arithmetic**

## Executive assessment

M1.1 is merged: x2/x3 have pinned structural verification, x3 includes its
proof-driven normalization repair, and x4 uses a cancellation-safe TwoSum +
`renormalize` baseline instead of the rejected FastTwoSum path.

M2 has now started with an independent higher-limb correctness oracle. The
current state is:

- x2 FMA structural proof: complete at the pinned verifier/toolchain;
- x3 FMA structural proof: complete after the final-normalization repair;
- x4 FMA: safe empirical baseline accepted; fixed-cost source-specific
  compression remains research;
- x5-x8 exact reference add/sub/mul/FMA: implemented under `Experimental` and
  cross-checked against 8192-bit packing;
- optimized/verified x5-x8 arithmetic networks: not implemented;
- formal proof of empirical constants `C_2`, `C_3`, and `C_4`: pending;
- MultiFloat-native DOT/SYRK/TRSM/GEMM backend: not implemented.

The package remains research software, but future x5-x8 optimization now has an
independent exact oracle instead of depending on another candidate arithmetic
network.

## M2 reference oracle

`MultiFloatArithmetic.Experimental` now provides:

- `reference_add`;
- `reference_sub`;
- `reference_mul`;
- `reference_fma`.

For Float32/Float64 expansions with `N in 5:8`, the scalar oracle converts each
normalized finite input to exact `Rational{BigInt}` dyadics, performs the
operation exactly, then packs the result once into an N-limb MultiFloat. Vector
reference operations are lane-wise scalar calls so future SIMD candidates do not
share implementation failure modes with their oracle.

The internal packing precision is not trusted by construction alone. CI compares
reference outputs to an independent 8192-bit BigFloat pack, checks normalization,
exact identities/cancellation, commutativity where required, Vec2 lane
semantics, explicit domain errors, and independence from ambient BigFloat
precision.

See `docs/REFERENCE_X5_X8.md` for the full oracle contract.

## x4 finding retained as a design constraint

The original x4 final compression began with
`fast_two_sum(b, a1)`. FPAN decomposition refuted the required magnitude
ordering under cancellation. In 10,000 seeded cases with `c` chosen as the x4
representation of `-x*y`:

- first-FastTwoSum ordering violations: **3,755 / 10,000**;
- old x4 non-normalized outputs: **5,258 / 10,000**;
- old output differed from the safe baseline: **5,258 / 10,000**;
- safe-baseline non-normalized outputs: **0 / 10,000**.

The M2/M3 rule is therefore explicit: higher-limb candidates may not inherit
FastTwoSum magnitude assumptions from lower-width code without proving the
source-specific ordering.

## Current x4 baseline

The correctness-first x4 path:

1. uses general `two_sum` instead of the two risky final-compression
   `fast_two_sum` operations;
2. finishes with `MultiFloats.renormalize`;
3. preserves the existing x4 product/accumulation network before final
   compression.

This gives up fixed-cost branch-free normalization for x4. Short
assumption-free fixed-tail A/B candidates did not prove the top two non-overlap
relations at p=24 or p=53.

## Performance snapshot

For the latest recorded safe-x4 generic hosted runner, Julia 1.10.11:

| Workload | Old x4 | Safe x4 | Relative result |
|---|---:|---:|---:|
| scalar, 20k | 0.062 ms | 0.345 ms | safe/old 5.593x |
| Vec2, 20k lanes | 0.155 ms | 0.211 ms | safe/old 1.364x |
| Vec4, 20k lanes | 0.080 ms | 0.109 ms | safe/old 1.352x |
| Vec8, 20k lanes | 0.064 ms | 0.068 ms | safe/old 1.056x |

Against upstream `x*y+c`, safe x4 SIMD measured about **1.59x (Vec2), 1.59x
(Vec4), and 1.31x (Vec8)** on that runner. Scalar x4 is not a performance
candidate.

The new M2 rational reference path is intentionally excluded from performance
comparison; it is an oracle, not a hot-loop implementation.

## Kernel decisions

| Component | Decision |
|---|---|
| x2 scalar `fma_fast` | marginal architecture-dependent opt-in candidate |
| x3 scalar `fma_fast` | structurally verified but slower; do not auto-select |
| x2/x3 SIMD `fma_fast` | continue architecture-specific optimization path |
| x4 scalar `fma_fast` | correctness baseline only; upstream is faster |
| x4 SIMD `fma_fast` | safe baseline retains useful speedup |
| x5-x8 `reference_*` | accepted only as exact Experimental oracle |
| optimized x5-x8 | next research stage; none accepted yet |
| quotient-digit division | rejected as default; retain under `Experimental` |
| expansion × one-limb multiply | experimental helper only |
| native linear algebra backend | not implemented |

## Acceptance gates

A higher-limb kernel may be promoted only when all applicable gates pass:

1. domain, rounding, underflow, normalization, and error metric are explicit;
2. exact identities and permanent adversarial cases agree with the M2 oracle;
3. normalized output is demonstrated and then proved where feasible;
4. every FastTwoSum magnitude precondition is proved or the operation is avoided;
5. discarded-tail/error constants have reviewable proof evidence;
6. scalar/SIMD lane semantics and required symmetries pass;
7. code generation is inspected only after correctness;
8. performance benefit is repeatable on the intended architecture/workload;
9. downstream solver A/B does not degrade residuals or certificates.

## Immediate next milestone

Finish M2 by merging and freezing the exact x5-x8 oracle. Then start M3 with a
**safe x5 addition baseline** that uses explicit normalization/discarded-tail
measurement and is differential-tested against `reference_add` before any
fixed-cost minimization. Extend to x6-x8 only after x5 has a defensible error
contract.
