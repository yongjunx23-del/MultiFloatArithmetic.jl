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
- Record first maximum safe-add constants: ~0.05283 (x5), ~0.02548 (x6),
  ~0.01099 (x7), and ~0.00542 (x8).
- Start M4 with `mul5_safe`, preserving every TwoProd product/residual and
  requiring pair-level exact reconstruction plus exact full-product tail
  accounting.
- Generalize multiplication to one `mul_safe` implementation for Float64x5-x8.
  Width N keeps all `2N^2` product/residual components, canonicalizes operand
  order, fully renormalizes, and truncates only after the exact product expansion
  is formed.
- Add x6-x8 pairwise TwoProd exactness, exact head+tail reconstruction,
  `reference_mul` equality, normalization, commutativity, dense/scaled/boundary
  tests, and separate empirical relative-error gates.
- Record first maximum `|err|/(u^N|x*y|)` values: ~0.04841 (x5), ~0.01837 (x6),
  ~0.00726 (x7), and ~0.00254 (x8); all measured diagnostics had zero
  oracle/normalization/commutativity failures.
- Record first safe-mul timing snapshots on Zen 3: 32.84 ms/500 (x5),
  5.059 ms/40 (x6), 5.710 ms/25 (x7), and 4.492 ms/12 (x8). Different case
  counts reflect rapidly growing `2N^2` exact-expansion cost; these are not
  performance targets.
