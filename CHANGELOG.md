# Changelog

## Unreleased

- Freeze the accepted top-level API around `fma_fast`.
- Move quotient-digit division and one-limb multiplication into the explicit
  `MultiFloatArithmetic.Experimental` namespace.
- Record the negative division A/B result instead of presenting the candidate as
  an optimization.
- Add the numerical contract, milestone status, benchmark interpretation, API
  boundary tests, binary32 smoke coverage, exact-identity checks, bounds-checked
  CI, and a macOS correctness job.
