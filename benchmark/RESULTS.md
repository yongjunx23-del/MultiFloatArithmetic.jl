# Recorded benchmark results

These are informational GitHub-hosted-runner snapshots, not portable hard gates.

## Safe Float64x5 multiplication — 2026-08-12 — Zen 3

Julia 1.10.11. `mul5_safe` is the M4 correctness baseline: all 25 limb products
are evaluated with TwoProd, both product/residual components are retained,
canonicalized, fully renormalized as 50 terms, and only then truncated to five
limbs.

| Corpus | Cases | Oracle mismatches | Normalization failures | Commutativity failures | Max `|err|/(u^5|xy|)` |
|---|---:|---:|---:|---:|---:|
| dense ordinary | 300 | 0 | 0 | 0 | 0.0484069 |
| scaled | 300 | 0 | 0 | 0 | 0.0417276 |
| powers of two | 25 | 0 | 0 | 0 | 0 |

Every dense/scaled case had a nonzero discarded tail while the five-limb result
still matched `reference_mul` bit-for-bit. The permanent unit corpus additionally
checks every individual TwoProd pair with exact `Rational{BigInt}` arithmetic.

First scalar timing:

- **33.557 ms for 500 Float64x5 products**.

This timing is intentionally poor and is not a regression target: sorting and
full 50-term renormalization exist to provide a defensible multiplication source
of truth before fixed-cost search. CI uses C=1 as an empirical relative-error
gate, not a proved bound.

The current safe-mul domain excludes nonfinite TwoProd components, nonzero
subnormal product/residual components, and pair products that underflow to zero.

## Safe Float64x5-x8 addition — 2026-08-12 — Zen 3

Julia 1.10.11. The implementation is the common M3 correctness baseline, not a
performance candidate.

| Width | Max observed `|err|/(u^N(|x|+|y|))` | Time / 5000 scalar adds |
|---|---:|---:|
| x5 | 0.0528336 | 0.454 ms |
| x6 | 0.0254805 | 0.545 ms |
| x7 | 0.0109878 | 0.710 ms |
| x8 | 0.00541727 | 1.438 ms |

Every measured width had zero oracle mismatches, normalization failures, and
commutativity failures. The near-cancellation diagnostic had zero discarded tail
for x5-x8. CI's C=1 bounds are empirical regression gates, not proofs.

## FMA context

Earlier snapshots established that scalar x3 and safe x4 custom FMA are not
performance wins, while SIMD x3 and safe x4 can retain architecture-dependent
benefits. The historical old x4 candidate is rejected because cancellation
violated a FastTwoSum precondition and produced non-normalized expansions.

The permanent 10k x4 cancellation diagnostic found 3,755 first-FastTwoSum
ordering violations and 5,258 old non-normalized outputs versus zero for the safe
baseline.

## Quotient-digit division

Specialized scalar quotient-digit division remains clearly noncompetitive;
vector results are architecture/type dependent. It remains Experimental only.

## Interpretation

- correctness evidence overrides wall-time wins;
- safe add5-add8 and mul5 are oracle-matching baselines, not optimized APIs;
- addition and multiplication keep operation/width-specific empirical tail
  measurements;
- fixed-cost minimization starts only after formal tail analysis;
- the current multiplication baseline is deliberately much slower than any
  future acceptable kernel;
- all performance claims require deployment-CPU and downstream solver A/B.
