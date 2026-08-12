# Project status

Status date: **2026-08-12**
Current milestone: **M4 baseline complete — safe Float64x5-x8 multiplication**

## Executive assessment

The project now has independent exact-oracle and correctness-first arithmetic
baselines through eight Float64 limbs:

- x2/x3 fixed-cost FMA networks with pinned structural verification;
- x3 proof-driven final-normalization repair;
- x4 cancellation-safe TwoSum + `renormalize` FMA baseline;
- exact-rational x5-x8 Experimental add/sub/mul/FMA oracles;
- one common no-FastTwoSum safe addition family for Float64x5-x8;
- one common no-FastTwoSum safe multiplication family for Float64x5-x8.

Formal discarded-tail proofs, fixed-cost high-limb networks, direct x5-x8 FMA,
and native linear algebra remain future work.

## M4 safe multiplication family

`Experimental.mul_safe` and wrappers `mul5_safe` ... `mul8_safe` use the same
algorithm for width N:

1. evaluate all N^2 limb pairs with `two_prod`;
2. retain every rounded product and residual, giving `2N^2` components;
3. reject nonfinite/subnormal/underflowing pair components under the current
   conservative domain;
4. canonicalize signed zero and sort by deterministic `(abs(value), value)` key;
5. fully renormalize the exact `2N^2`-term expansion;
6. retain the N-limb head and renormalize it.

No FastTwoSum is used.

The permanent x5-x8 corpora check every accepted TwoProd pair with exact
`Rational{BigInt}` arithmetic and require:

```text
value(result) + value(discarded normalized tail) = value(x) * value(y)
```

They also require canonical term equality under operand swap, normalized full
and returned expansions, bitwise `reference_mul` equality, bitwise
commutativity, exact identities/sign cases, dense/scaled operands, and
power/boundary products.

All M4 width-specific unit gates are green on Linux Julia 1.10/current and macOS
current.

### First width-specific relative-error measurements

Zen 3, Julia 1.10.11:

| Width | Max observed `|err|/(u^N|xy|)` | Informational scalar timing |
|---|---:|---:|
| x5 | 0.0484069 | 32.84 ms / 500 |
| x6 | 0.0183703 | 5.059 ms / 40 |
| x7 | 0.00725543 | 5.710 ms / 25 |
| x8 | 0.00253803 | 4.492 ms / 12 |

All dense/scaled diagnostics had zero oracle, normalization, and commutativity
failures. Power-of-two boundary cases were exact. CI keeps a width-specific
`C=1` empirical gate; it is not a theorem.

The timings are intentionally poor and use different case counts. The safe family
allocates, sorts up to 128 components, and fully renormalizes them; it is a
correctness source of truth, not a performance API.

## Current multiplication domain

The baseline accepts normalized finite Float64 N=5:8 values but deliberately
rejects pair computations with:

- nonfinite TwoProd product/residual components;
- nonzero subnormal products;
- nonzero subnormal TwoProd residuals;
- nonzero input pairs whose rounded product underflows to zero.

This remains conservative until gradual-underflow semantics and exact TwoProd
behavior are analyzed there.

## M3 safe addition family

`Experimental.add_safe` and `add5_safe` ... `add8_safe` use N general TwoSums,
full exact 2N-term renormalization, N-limb truncation, and head renormalization.
No FastTwoSum is used.

First max observed `|err|/(u^N(|x|+|y|))` values were 0.0528336 (x5),
0.0254805 (x6), 0.0109878 (x7), and 0.00541727 (x8), all below the empirical
C=1 gates.

## Current decisions

| Component | Decision |
|---|---|
| x2/x3 fixed-cost FMA | structurally verified; scalar speedups are architecture-sensitive |
| x4 FMA | cancellation-safe correctness baseline; SIMD still promising |
| x5-x8 `reference_*` | exact Experimental oracle |
| x5-x8 `add_safe` | accepted M3 Experimental correctness family |
| x5-x8 `mul_safe` | accepted M4 Experimental correctness family |
| fixed-cost high-limb add/mul | blocked on formal tail analysis |
| direct x5-x8 FMA/submul | M5 next |
| reciprocal/div/sqrt | M6 |
| native DOT/SYRK/TRSM/GEMM | M7 |

## Immediate next milestone

Start **M5 with a correctness-first Float64x5 direct FMA baseline**. Build an
over-complete exact `x*y+c` expansion using the accepted multiplication
components plus exact addend components, canonicalize and fully renormalize it,
then require exact head+tail accounting and bitwise agreement with
`reference_fma` before any fused product-network minimization.
