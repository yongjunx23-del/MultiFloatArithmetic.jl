# Verified Float64x5-Float64x8 arithmetic research plan

Research-only until formal verification and numerical gates pass.

## Acceptance contract

A future performance kernel must have normalized output, explicit discarded-tail
bounds, required symmetries, adversarial exact-oracle agreement, proven
FastTwoSum assumptions, and reproducible IEEE round-to-nearest behavior before
branch-free minimization is considered.

The x4 cancellation defect permanently sets the search order: **correctness and
normalization before gate count**.

## M2 exact oracle — complete baseline

`Experimental.reference_add/sub/mul/fma` provides the independent x5-x8 rejection
standard via exact `Rational{BigInt}` arithmetic and one final pack. CI
cross-checks the pack against independent 8192-bit BigFloat results.

## M3 addition — safe family complete baseline

`Experimental.add_safe` covers Float64 x5-x8 with N same-index general TwoSums,
full exact 2N-term renormalization, N-limb truncation, and head renormalization.
No FastTwoSum is used.

First maximum observed `|err|/(u^N(|x|+|y|))` values are ~0.0528336 (x5),
~0.0254805 (x6), ~0.0109878 (x7), and ~0.00541727 (x8). C=1 remains an
empirical regression gate, not a proof.

## M4 multiplication — safe family complete baseline

`Experimental.mul_safe` covers Float64 x5-x8. For width N:

1. compute every N×N limb pair with TwoProd;
2. preserve all rounded products and residuals (`2N^2` components);
3. reject nonfinite/subnormal/underflowing pair components under the current
   conservative domain;
4. canonicalize signed zero and sort components so operand swap yields an
   identical accumulation sequence;
5. fully renormalize the exact component tuple;
6. retain the N-limb head and renormalize it.

No FastTwoSum is used.

Permanent gates require exact pairwise TwoProd decomposition, exact returned-head
+ discarded-tail reconstruction, bitwise `reference_mul` agreement,
normalization, commutativity, exact identities/signs, and dense/scaled/boundary
corpora.

First maximum observed `|err|/(u^N|xy|)` values are:

- x5: ~0.0484069
- x6: ~0.0183703
- x7: ~0.00725543
- x8: ~0.00253803

All measured width-specific diagnostics have zero oracle/normalization/
commutativity failures. C=1 remains an empirical gate.

### Add/mul proof work before minimization

Derive reviewable discarded-tail bounds for the common safe addition and
multiplication constructions. Do not search fixed-cost add/mul networks merely
because the empirical constants are small. Every component reordering or
FastTwoSum substitution must preserve exact-oracle and symmetry gates and have a
source-specific proof argument.

## M5 direct FMA — next implementation milestone

Start with Float64x5 only. The first candidate should form an over-complete exact
`x*y+c` expansion from:

- all 25 TwoProd product/residual pairs used by M4;
- all five addend limbs `c[i]` as exact components;
- deterministic canonical ordering before full renormalization.

Then retain the five-limb head only after the complete 55-term expansion has
been normalized.

Acceptance sequence:

`fma5_safe -> exact component/head-tail/oracle/symmetry gates -> tail study ->
safe fma6/fma7/fma8 -> formal bound -> fixed-cost direct-FMA search`.

`fma5_safe` must compare directly with `reference_fma`, not with `mul_safe` plus
`add_safe`, so the oracle remains implementation-independent.

Do not use the existing x2-x4 `fma_fast` networks as structural templates for
x5; their fixed-cost assumptions belong to a different proof stage.

## M6-M7 direction

After safe/verified direct FMA exists:

- M6: reciprocal/division/sqrt via 1->2->4->8 precision doubling, using accepted
  add/mul/FMA primitives and stronger result-relative paths where needed;
- M7: MultiFloatVec-native DOT/SYRK/TRSM/GEMM, then TRSM/Cholesky and downstream
  solver A/B on residuals, certificates, time, and memory.

Every discovered counterexample becomes a permanent regression test.
