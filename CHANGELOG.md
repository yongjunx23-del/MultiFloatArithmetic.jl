# Changelog

## Unreleased

- Freeze the top-level research API around `fma_fast` and move rejected
  quotient-digit division / one-limb multiplication into
  `MultiFloatArithmetic.Experimental`.
- Add the numerical contract, milestone status, architecture-dependent benchmark
  records, API-boundary tests, binary32 smoke coverage, exact-identity checks,
  bounds-checked Linux/macOS CI, and code-generation diagnostics.
- Add pinned FPANVerifier/Z3 structural-proof assets for the fixed-cost x2/x3
  source networks.
- Repair x3 final normalization with the minimal proved TwoSum/FastTwoSum pair.
- Identify a cancellation-invalid first FastTwoSum in the original x4 end
  network and reproduce the problem on concrete Float64 inputs.
- Add the correctness-first x4 baseline: general TwoSum final compression plus
  `MultiFloats.renormalize`.
- Add an old-versus-safe x4 scalar/SIMD benchmark and exact-cancellation
  diagnostic. In the seeded 10k cancellation corpus the old path produced 5,258
  non-normalized results while the safe baseline produced none.
- Start M2 with exact-rational x5-x8 Experimental reference add/sub/mul/FMA for
  Float32/Float64, including lane-wise `MultiFloatVec` wrappers.
- Cross-check the x5-x8 reference packing precision against independent 8192-bit
  BigFloat packing and add explicit domain/normalization/ambient-precision tests.
