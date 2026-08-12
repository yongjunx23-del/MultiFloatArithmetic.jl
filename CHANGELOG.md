# Changelog

## Unreleased

- Freeze the top-level research API around `fma_fast`.
- Move quotient-digit division and one-limb multiplication into the explicit
  `MultiFloatArithmetic.Experimental` namespace.
- Record the mixed, architecture-dependent division A/B evidence instead of
  presenting the candidate as a general optimization.
- Add the numerical contract, milestone status, two-architecture benchmark
  record, API-boundary tests, binary32 smoke coverage, exact-identity checks,
  bounds-checked CI, and a macOS correctness job.
- Add pinned FPANVerifier/Z3 structural-proof assets for the x2-x4 source
  networks.
- Repair the x3 final normalization with the minimal proved TwoSum/FastTwoSum
  pair after the original leading non-overlap obligation was refuted.
- Add an old-versus-repaired x3 scalar/SIMD benchmark so the proof repair's fixed
  cost remains visible.
