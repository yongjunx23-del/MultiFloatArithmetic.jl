# Changelog

## Unreleased

- Freeze the accepted top-level API around `fma_fast`.
- Move quotient-digit division and one-limb multiplication into the explicit
  `MultiFloatArithmetic.Experimental` namespace.
- Record the mixed, architecture-dependent division A/B evidence instead of
  presenting the candidate as a general optimization.
- Add the numerical contract, milestone status, two-architecture benchmark
  record, API-boundary tests, binary32 smoke coverage, exact-identity checks,
  bounds-checked CI, and a macOS correctness job.
