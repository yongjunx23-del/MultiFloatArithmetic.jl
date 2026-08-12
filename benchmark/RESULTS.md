# Recorded benchmark results

These are informational GitHub-hosted-runner snapshots, not portable hard gates.

## Safe Float64x5-x8 multiplication — 2026-08-12 — Zen 3

Julia 1.10.11. `mul_safe` is the M4 correctness family, not a performance
candidate. Width N retains and fully normalizes all `2N^2` TwoProd product and
residual components before N-limb truncation.

| Width | Max observed `|err|/(u^N|xy|)` | Informational timing |
|---|---:|---:|
| x5 | 0.0484069 | 32.84 ms / 500 |
| x6 | 0.0183703 | 5.059 ms / 40 |
| x7 | 0.00725543 | 5.710 ms / 25 |
| x8 | 0.00253803 | 4.492 ms / 12 |

For x6, dense/scaled diagnostics used 80+80 cases; for x7, 60+60; for x8,
40+40. Every dense/scaled case had a nonzero discarded tail and every measured
width had zero oracle mismatches, normalization failures, and commutativity
failures. Power-of-two boundary cases were exact.

The permanent unit corpus also checks each individual TwoProd pair with exact
`Rational{BigInt}` reconstruction and verifies exact returned-head + discarded-
tail reconstruction of the full product.

CI uses a width-specific C=1 empirical relative-error gate. No value above is a
formal worst-case bound.

The timing columns use different case counts and must not be compared as raw
throughputs without normalization. They only show that the over-complete safe
family is intentionally expensive, especially as x8 reaches 128 exact product
components.

## Safe Float64x5-x8 addition — 2026-08-12 — Zen 3

| Width | Max observed `|err|/(u^N(|x|+|y|))` | Time / 5000 scalar adds |
|---|---:|---:|
| x5 | 0.0528336 | 0.453 ms |
| x6 | 0.0254805 | 0.544 ms |
| x7 | 0.0109878 | 0.709 ms |
| x8 | 0.00541727 | 1.442 ms |

Every measured width had zero oracle mismatches, normalization failures, and
commutativity failures. Near-cancellation diagnostics had zero discarded tail.

## FMA context

Scalar x3 and safe x4 custom FMA are not performance wins on the recorded
runners, while SIMD x3 and safe x4 retain architecture-dependent benefits. The
historical old x4 candidate remains rejected because cancellation violated a
FastTwoSum precondition and produced non-normalized expansions.

## Quotient-digit division

Specialized scalar quotient-digit division remains noncompetitive; vector
results are architecture/type dependent. It remains Experimental only.

## Interpretation

- correctness evidence overrides wall-time wins;
- safe x5-x8 add/mul are oracle-matching baselines, not optimized APIs;
- each operation and width keeps its own empirical tail measurement;
- formal tail analysis precedes fixed-cost minimization;
- all future performance claims require deployment-CPU and downstream solver A/B.
