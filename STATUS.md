# Project status

Status date: **2026-08-12**
Current milestone: **M3 baseline complete — safe Float64x5-x8 addition**

## Executive assessment

The project now has a correctness-first arithmetic ladder:

- x2/x3 FMA fixed-cost networks with pinned structural verification;
- x3 proof-driven final-normalization repair;
- x4 cancellation-safe TwoSum + `renormalize` baseline;
- exact-rational x5-x8 Experimental add/sub/mul/FMA oracles;
- one common no-FastTwoSum safe addition construction accepted as the M3
  Experimental baseline for Float64x5, x6, x7, and x8.

Optimized higher-limb arithmetic, formal discarded-tail proofs, and native linear
algebra remain future work.

## M3 safe addition family

`Experimental.add_safe` and wrappers `add5_safe` ... `add8_safe` use the same
algorithm at every width N:

1. N same-index general TwoSums;
2. full renormalization of the exact 2N-term expansion;
3. retain the N-limb head;
4. renormalize the head.

CI requires exact head+tail accounting, normalized output, bitwise
commutativity, and bitwise equality with the M2 `reference_add` oracle under
ordinary, wide-exponent, exact-identity, deep-cancellation, power/carry-boundary,
and unbalanced-operand corpora.

All x5-x8 safe-add unit gates are green on Linux Julia 1.10/current and macOS
current.

### First width-specific tail measurements

| Width | Max observed `|err|/(u^N(|x|+|y|))` | First timing / 5000 |
|---|---:|---:|
| x5 | 0.0528336 | 0.454 ms |
| x6 | 0.0254805 | 0.545 ms |
| x7 | 0.0109878 | 0.710 ms |
| x8 | 0.00541727 | 1.438 ms |

Near-cancellation diagnostics had zero discarded tail at all four widths in the
first accepted corpus. CI uses C=1 as a conservative **empirical regression
gate**, not a theorem. See `docs/ADD_SAFE_X5_X8.md`.

## M2 exact oracle

`Experimental.reference_add`, `reference_sub`, `reference_mul`, and
`reference_fma` support Float32/Float64 N=5:8. Scalar operations use exact
`Rational{BigInt}` dyadics and pack once; CI independently cross-checks packing
against 8192-bit BigFloat. Vector reference methods are deliberately lane-wise.

## x4 design constraint

The rejected old x4 `fast_two_sum(b,a1)` assumption failed under cancellation.
In 10,000 seeded exact-cancellation-style cases, 3,755 violated its ordering and
5,258 old outputs were non-normalized; the safe x4 baseline produced none.
Higher-limb work therefore cannot inherit FastTwoSum assumptions without proof.

## Current decisions

| Component | Decision |
|---|---|
| x2 scalar FMA | marginal architecture-dependent candidate |
| x3 scalar FMA | structurally verified but slower; no auto-selection |
| x2/x3 SIMD FMA | continue architecture-specific path |
| x4 scalar FMA | correctness baseline only |
| x4 SIMD FMA | safe baseline retains useful speedup |
| x5-x8 `reference_*` | exact Experimental correctness oracle |
| x5-x8 `add_safe` | accepted M3 Experimental correctness baseline |
| fixed-cost add5-add8 | not started; formal tail analysis first |
| safe multiplication x5 | next implementation milestone |
| quotient-digit division | rejected as default; Experimental only |
| native linear algebra | not implemented |

## Acceptance gates for higher-limb kernels

1. explicit domain/rounding/normalization/error contract;
2. permanent adversarial cases agree with the exact M2 oracle;
3. exact discarded-tail accounting where available;
4. normalized outputs and required symmetries;
5. every FastTwoSum precondition proved or FastTwoSum avoided;
6. empirical constants measured, then replaced or backed by formal bounds;
7. SIMD lane semantics where applicable;
8. codegen/performance only after correctness;
9. downstream solver A/B must preserve residuals and certificates.

## Immediate next milestone

M3's safe family is ready to freeze. In parallel with deriving a formal
higher-limb addition tail bound, start **M4 with a correctness-first Float64x5
multiplication baseline**: over-complete product expansion, general error-free
transforms/renormalization, exact discarded-tail accounting, and differential
comparison against `reference_mul` before any FastTwoSum or performance search.
