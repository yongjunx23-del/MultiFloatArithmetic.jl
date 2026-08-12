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
- Start M5 with `Experimental.fma5_safe`, a direct Float64x5 55-component
  `x*y+c` correctness baseline built from all 25 TwoProd product/residual pairs
  plus the five addend limbs before any five-limb truncation.
- Add exact pairwise product checks, exact 55-term FMA head+tail reconstruction,
  `reference_fma` equality, x/y symmetry, dense/scaled identities, powers of two,
  and destructive-cancellation regression cases for direct x5 FMA.
- Add an empirical direct-FMA C=1 gate for `|err|/(u^5(|xy|+|c|))`; the first
  diagnostics measured maxima ~0.04568 ordinary and ~0.04928 scaled with zero
  direct oracle/normalization/symmetry failures.
- Compare direct FMA to the rounded composition `add_safe(mul_safe(x,y),c)`.
  The composition disagreed with `reference_fma` in 45/200 ordinary, 41/200
  scaled, and 150/150 destructive-cancellation cases, while the direct path had
  zero oracle mismatches.
- Record the first x5 timing A/B: 7.707 ms / 120 direct safe FMA operations vs
  6.231 ms / 120 safe mul-then-add operations on the cited Zen 3 runner. The
  direct correctness baseline is currently slower but avoids the intermediate
  rounding error that dominates cancellation-sensitive use.
