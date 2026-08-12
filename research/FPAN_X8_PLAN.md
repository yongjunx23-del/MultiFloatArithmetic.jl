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

A single `Experimental.add_safe` construction covers Float64 x5-x8:

1. N same-index general TwoSums;
2. full renormalization of the exact 2N terms;
3. retain N leading terms and renormalize the head.

Width-specific wrappers `add5_safe` ... `add8_safe` share this implementation.
No FastTwoSum is used.

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

## M4 multiplication — Float64x5 safe baseline accepted

`Experimental.mul5_safe` is deliberately over-complete:

1. compute every 5x5 limb pair with TwoProd;
2. preserve all 25 rounded products and all 25 residuals;
3. canonicalize signed zeros and sort the 50 finite components so operand swap
   produces an identical term sequence;
4. fully renormalize the exact 50-term expansion;
5. retain five leading limbs and renormalize the head.

No FastTwoSum is used.

Permanent gates require each individual TwoProd to reconstruct its exact dyadic
pair product, exact returned-head + discarded-tail accounting, bitwise
`reference_mul` agreement, normalized output/full expansion, and bitwise
commutativity.

The first diagnostics found zero oracle/normalization/commutativity failures and
max `|err|/(u^5|xy|) ≈ 0.0484069`. CI uses C=1 only as an empirical gate.
The current baseline rejects nonfinite or subnormal TwoProd components and pair
products that underflow to zero until underflow semantics are proved.

The first scalar timing was 33.557 ms for 500 products on Zen 3. This cost is
acceptable only because M4 is still establishing a correctness source of truth.

### M4 next

Generalize the same over-complete product construction to safe Float64x6, x7,
and x8. Each width must independently establish:

- exact primitive TwoProd decomposition in its accepted domain;
- exact full-product head+tail accounting;
- bitwise `reference_mul` equality;
- commutativity and normalization;
- separate empirical relative-error constants;
- explicit overflow/subnormal/underflow behavior.

Acceptance sequence:

`mul5_safe -> safe mul6/mul7/mul8 -> formal multiplication tail bound -> fixed-cost search`.

Do not copy lower-width FastTwoSum assumptions into M4. Do not optimize away the
canonicalization until a proposed accumulation ordering has a commutativity and
error proof.

## M5-M7 direction

After safe/verified add and mul exist:

- M5: direct FMA/submul x5-x8, then DOT/SYRK/GEMM hot paths;
- M6: reciprocal/division/sqrt via 1->2->4->8 precision doubling;
- M7: MultiFloatVec-native DOT/SYRK/TRSM/GEMM and downstream solver A/B.

Every discovered counterexample becomes a permanent regression test.
