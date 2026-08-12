# Changelog

## Unreleased

- Freeze the top-level research API around `fma_fast`; keep rejected quotient
  division and one-limb multiplication under `Experimental`.
- Add numerical contracts, cross-platform bounds-checked CI, benchmark records,
  API/identity/cancellation tests, and code-generation diagnostics.
- Add pinned FPANVerifier/Z3 structural proof assets for fixed-cost x2/x3 FMA.
- Repair x3 final normalization with the minimal proved TwoSum/FastTwoSum pair.
- Identify and reproduce the cancellation-invalid FastTwoSum in the old x4 path;
  replace it with the correctness-first TwoSum + `renormalize` baseline.
- Add exact-rational x5-x8 Experimental reference add/sub/mul/FMA and independent
  8192-bit packing checks.
- Complete M3 with one no-FastTwoSum `add_safe` implementation for Float64x5-x8,
  exact tail/oracle/commutativity gates, and width-specific empirical C=1 gates.
- Complete M4 with one no-FastTwoSum `mul_safe` implementation for Float64x5-x8,
  exact pairwise TwoProd and full-product tail accounting, oracle equality,
  commutativity, conservative underflow semantics, and width-specific gates.
- Complete M5 with one direct `fma_safe` family for Float64x5-x8. Width N combines
  all `2N^2` exact TwoProd product/residual components with all N addend limbs
  before N-limb truncation; full component counts are 55/78/105/136 for x5-x8.
- Require exact pairwise product checks, exact direct-FMA head+tail reconstruction,
  `reference_fma` equality, x/y symmetry, dense/scaled identities, powers of two,
  and destructive-cancellation regression cases at every x5-x8 width.
- Record first maximum direct-FMA constants in `u^N(|xy|+|c|)` units:
  ~0.04928 (x5), ~0.02276 (x6), ~0.00452 (x7), and ~0.000826 (x8), with zero
  direct oracle/normalization/symmetry failures.
- Confirm that rounded `add_safe(mul_safe(x,y),c)` is not an exact substitute for
  direct FMA: it disagreed with `reference_fma` in all destructive-cancellation
  cases at x5 (150/150), x6 (48/48), x7 (32/32), and x8 (24/24).
- Record first safe direct-vs-composed timing ratios of 0.811x, 0.837x, 0.856x,
  and 0.882x (composition/direct) for x5-x8. The direct correctness path remains
  slower and is not yet a hot kernel.
- Add a permanent BigFloat/MultiFloat precision-throughput benchmark for the
  requested BigFloat128/x2, BigFloat256/x4, and BigFloat512/x8 pairings. Original
  upstream x2/x4 arithmetic strongly outperformed BigFloat128/256 in the first
  hosted snapshot, while the intentionally over-complete safe x8 multiplication
  was roughly 3,600x slower than BigFloat512, making x8 network reduction a
  priority before native linear algebra.
