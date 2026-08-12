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

A single `Experimental.add_safe` construction now covers Float64 x5-x8:

1. N same-index general TwoSums;
2. full renormalization of the exact 2N terms;
3. retain N leading terms and renormalize the head.

Width-specific wrappers `add5_safe` ... `add8_safe` share this implementation.
No FastTwoSum is used.

Permanent width-specific corpora require exact head+tail accounting, bitwise
`reference_add` agreement, normalization, commutativity, exact identities, deep
cancellation, power/carry boundaries, and unbalanced operands.

First maximum observed `|err|/(u^N(|x|+|y|))` values:

- x5: ~0.0528336
- x6: ~0.0254805
- x7: ~0.0109878
- x8: ~0.00541727

All current empirical CI gates use C=1 and are explicitly non-theorem.

### Addition proof/minimization next

Derive a reviewable discarded-tail bound for the common 2N-term construction.
Only after that proof work should fixed-cost add5-add8 networks be searched.
Every proposed FastTwoSum must prove its source-specific ordering.

## M4 multiplication — next implementation milestone

Start with Float64x5 only. The first candidate should be deliberately
over-complete: retain enough exact product components/error terms to account for
the full product before final N-limb truncation, fully normalize, compare against
`reference_mul`, and measure the discarded normalized tail.

Acceptance sequence:

`mul5_safe -> exact-tail/oracle/cancellation/commutativity gates -> tail-bound
study -> safe mul6/mul7/mul8 -> formal bound -> fixed-cost search`.

Do not copy lower-width FastTwoSum assumptions into M4.

## M5-M7 direction

After safe/verified add and mul exist:

- M5: direct FMA/submul x5-x8, then DOT/SYRK/GEMM hot paths;
- M6: reciprocal/division/sqrt via 1->2->4->8 precision doubling;
- M7: MultiFloatVec-native DOT/SYRK/TRSM/GEMM and downstream solver A/B.

Every discovered counterexample becomes a permanent regression test.
