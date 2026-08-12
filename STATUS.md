# Project status

Status date: **2026-08-12**
Current milestone: **M4 — safe Float64x5 multiplication baseline accepted**

## Executive assessment

The project now has a correctness-first arithmetic ladder:

- x2/x3 FMA fixed-cost networks with pinned structural verification;
- x3 proof-driven final-normalization repair;
- x4 cancellation-safe TwoSum + `renormalize` baseline;
- exact-rational x5-x8 Experimental add/sub/mul/FMA oracles;
- one common no-FastTwoSum safe addition construction accepted for Float64x5-x8;
- a correctness-first Float64x5 multiplication baseline that preserves every
  TwoProd component before final five-limb truncation.

Optimized higher-limb arithmetic, formal discarded-tail proofs, and native linear
algebra remain future work.

## M4 Float64x5 safe multiplication

`Experimental.mul5_safe` is deliberately over-complete:

1. evaluate all 25 limb-pair products with `two_prod`;
2. retain both rounded product and residual for every pair, giving 50 terms;
3. canonicalize signed zeros and sort the finite terms by deterministic
   magnitude/value key so operand swap produces the same term sequence;
4. fully renormalize the 50-term expansion;
5. retain the leading five limbs and renormalize the head.

No FastTwoSum is used.

The permanent corpus checks each individual TwoProd decomposition with exact
`Rational{BigInt}` arithmetic, then requires:

```text
value(result) + value(discarded 45 normalized limbs) = value(x) * value(y)
```

It also requires bitwise equality with `reference_mul`, normalized full/output
expansions, canonical-term equality under operand swap, and bitwise
commutativity. Dense multi-limb, scaled, sign, zero/one, power-of-two, and
near-boundary cases are covered.

The current M4 domain is intentionally conservative around gradual underflow:
nonfinite TwoProd components, nonzero subnormal product/residual components, and
pair-product underflow to zero are rejected until a dedicated underflow proof is
available.

### First multiplication diagnostic

Zen 3 hosted runner, Julia 1.10.11:

| Corpus | Cases | Oracle / norm / comm failures | Max `|err|/(u^5|xy|)` |
|---|---:|---:|---:|
| dense ordinary | 300 | 0 / 0 / 0 | 0.0484069 |
| scaled | 300 | 0 / 0 / 0 | 0.0417276 |
| powers of two | 25 | 0 / 0 / 0 | 0 |

Every nontrivial dense/scaled case had a nonzero discarded tail while still
matching the exact five-limb oracle bit-for-bit. The first scalar timing was
**33.557 ms for 500 products**. This cost is expected: sorting and full 50-term
renormalization make this a correctness baseline, not a performance candidate.

CI uses `C = 1` as a conservative empirical relative-error gate. It is not a
formal theorem. See `docs/MUL5_SAFE.md`.

## M3 safe addition family

`Experimental.add_safe` and wrappers `add5_safe` ... `add8_safe` use N general
TwoSums, full renormalization of the exact 2N terms, N-limb truncation, and head
renormalization. No FastTwoSum is used.

First maximum observed `|err|/(u^N(|x|+|y|))` values were:

| Width | Constant | Time / 5000 scalar adds |
|---|---:|---:|
| x5 | 0.0528336 | 0.454 ms |
| x6 | 0.0254805 | 0.545 ms |
| x7 | 0.0109878 | 0.710 ms |
| x8 | 0.00541727 | 1.438 ms |

All safe-add width-specific corpora are green on Linux Julia 1.10/current and
macOS current. CI's C=1 gates remain empirical, not proved bounds.

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
| `mul5_safe` | accepted M4 Experimental correctness baseline; deliberately slow |
| safe multiplication x6-x8 | next M4 step |
| fixed-cost higher-limb add/mul | formal tail analysis first |
| quotient-digit division | rejected as default; Experimental only |
| native linear algebra | not implemented |

## Acceptance gates for higher-limb kernels

1. explicit domain/rounding/normalization/error contract;
2. permanent adversarial cases agree with the exact M2 oracle;
3. exact primitive decomposition and discarded-tail accounting where available;
4. normalized outputs and required symmetries;
5. every FastTwoSum precondition proved or FastTwoSum avoided;
6. empirical constants measured, then replaced or backed by formal bounds;
7. SIMD lane semantics where applicable;
8. codegen/performance only after correctness;
9. downstream solver A/B must preserve residuals and certificates.

## Immediate next milestone

Generalize the M4 multiplication construction to Float64x6, x7, and x8 with a
separate oracle/error gate at each width. Keep the same rule as M3: first build
an over-complete exact-product baseline, then derive a formal discarded-tail
bound, and only then search fixed-cost multiplication networks or SIMD hot paths.
