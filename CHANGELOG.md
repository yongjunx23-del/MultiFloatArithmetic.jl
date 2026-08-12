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
- Add the M3 Float64x5 safe addition baseline with exact tail accounting and an
  empirical C=1 regression gate.
- Generalize the same no-FastTwoSum safe addition construction through x8 using
  one `add_safe` implementation plus width-specific add5-add8 wrappers.
- Add cross-platform x6-x8 oracle/normalization/commutativity/cancellation/
  power-boundary tests and independent width-specific tail diagnostics.
- Record first maximum observed addition constants: ~0.05283 (x5), ~0.02548
  (x6), ~0.01099 (x7), and ~0.00542 (x8), all below the empirical C=1 gate.
- Record first safe-add scalar timings for 5,000 operations: 0.454 ms (x5),
  0.545 ms (x6), 0.710 ms (x7), 1.438 ms (x8) on the cited Zen 3 runner.
