# Verified Float64x5-Float64x8 arithmetic research plan

Research-only until formal verification and numerical gates pass.

## Core rule

The x4 cancellation defect permanently sets the search order: correctness,
normalization, exact-oracle agreement, and explicit error accounting before
branch-free/fixed-cost minimization. Every FastTwoSum must have a proved
source-specific magnitude precondition.

## M2 exact oracle — complete baseline

`Experimental.reference_add/sub/mul/fma` provides the independent x5-x8 rejection
standard via exact `Rational{BigInt}` arithmetic and one final pack.

## M3 addition — safe family complete baseline

`Experimental.add_safe` covers Float64 x5-x8 with exact 2N-term accumulation and
full renormalization before N-limb truncation. First empirical constants in
`u^N(|x|+|y|)` units are ~0.05283, 0.02548, 0.01099, and 0.00542 for x5-x8.

## M4 multiplication — safe family complete baseline

`Experimental.mul_safe` covers Float64 x5-x8 by preserving all `2N^2` TwoProd
components, canonicalizing their order, fully renormalizing, then truncating.
First empirical relative constants are ~0.04841, 0.01837, 0.00726, and 0.00254.

## M5 direct FMA — Float64x5 safe baseline accepted

`Experimental.fma5_safe` deliberately does not call the rounded N-limb
multiplication result. It combines:

- all 50 exact x5 product/residual components from the 25 TwoProd pairs;
- all five exact limbs of `c`;
- deterministic canonical ordering;
- one full renormalization of the resulting 55 components;
- five-limb truncation only after the exact `x*y+c` expansion exists.

Permanent gates require exact TwoProd reconstruction, exact 55-term head+tail
accounting, normalized full/output expansions, `reference_fma` equality, x/y
symmetry, identities, dense/scaled cases, powers of two, and deep destructive
cancellation.

First direct-FMA empirical max `|err|/(u^5(|xy|+|c|))` values were ~0.04568
ordinary and ~0.04928 scaled, with zero direct oracle/normalization/symmetry
failures.

### Why direct FMA is now mandatory research, not optional syntax

The rounded composition `add_safe(mul_safe(x,y),c)` disagreed with the exact
`reference_fma` oracle in:

- 45 / 200 ordinary cases;
- 41 / 200 scaled cases;
- 150 / 150 destructive-cancellation cases.

The direct path had zero oracle mismatches. Therefore separately correct add/mul
baselines are not a substitute for direct FMA in cancellation-sensitive residual,
refinement, or certificate work.

The initial direct correctness path is ~24% slower than the safe composition on
the first Zen 3 timing sample; this is acceptable at the baseline stage.

### M5 next

Generalize to one safe direct-FMA family for Float64 x5-x8. Width N should combine
all `2N^2` product/residual components with N addend limbs, producing:

- x5: 55 components;
- x6: 78 components;
- x7: 105 components;
- x8: 136 components.

For x6-x8 independently require:

- exact pairwise TwoProd reconstruction;
- exact direct FMA head+tail accounting;
- bitwise `reference_fma` equality;
- x/y symmetry and normalized output;
- destructive-cancellation coverage;
- separate empirical `u^N(|xy|+|c|)` constants;
- direct-vs-rounded-composition mismatch telemetry;
- explicit inherited underflow/overflow domain.

Acceptance sequence:

`fma5_safe -> safe fma6/fma7/fma8 -> formal direct-FMA tail analysis -> fixed-cost fused search`.

Do not use the existing x2-x4 `fma_fast` networks as templates for x5-x8 fixed-
cost structure before the safe direct family is complete.

## M6-M7 direction

- M6: reciprocal/division/sqrt via precision doubling, with direct FMA/submul for
  residual corrections where it materially improves error propagation.
- M7: MultiFloatVec-native DOT/SYRK/GEMM, then TRSM/Cholesky and downstream solver
  A/B on iterations, residuals, certificates, time, and memory.

Every discovered counterexample becomes a permanent regression test.
