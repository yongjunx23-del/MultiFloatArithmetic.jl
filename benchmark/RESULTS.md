# Recorded benchmark results

These are informational GitHub-hosted-runner snapshots, not portable hard gates.

## Safe Float64x5 direct FMA — 2026-08-12 — Zen 3

Julia 1.10.11. `fma5_safe` is the first M5 correctness baseline. It forms and
fully renormalizes one 55-component exact `x*y+c` expansion before five-limb
truncation.

| Corpus | Cases | Direct oracle failures | Mul-then-add oracle failures | Max `|err|/(u^5(|xy|+|c|))` |
|---|---:|---:|---:|---:|
| dense ordinary | 200 | 0 | 45 | 0.0456757 |
| scaled | 200 | 0 | 41 | 0.0492768 |
| destructive cancellation | 150 | 0 | 150 | ~1.01e-67 |

The direct path also had zero normalization and x/y-symmetry failures. Its
maximum result-relative errors in these seeded corpora were of order `1e-81`.
The cancellation row is particularly important: the five-limb rounded
composition `add_safe(mul_safe(x,y),c)` disagreed with the exact direct-FMA
oracle in **every tested cancellation case**, despite both component operations
having their own exact-oracle baselines.

Timing for 120 cases:

- direct 55-component safe FMA: **7.707 ms**;
- safe `mul_safe` followed by `add_safe`: **6.231 ms**;
- composition/direct = **0.809x**.

Thus the first direct baseline is about 24% slower than the already-slow safe
composition. M5 currently values it for avoiding intermediate-rounding loss;
performance optimization begins only after the direct x5-x8 family and formal
error analysis exist.

CI uses C=1 for the operand-relative direct-FMA metric as an empirical gate, not
a theorem.

## Safe Float64x5-x8 multiplication — 2026-08-12 — Zen 3

| Width | Max observed `|err|/(u^N|xy|)` | Informational timing |
|---|---:|---:|
| x5 | 0.0484069 | 25.734 ms / 500 |
| x6 | 0.0183703 | 4.062 ms / 40 |
| x7 | 0.00725543 | 4.578 ms / 25 |
| x8 | 0.00253803 | 3.626 ms / 12 |

Every measured width had zero oracle, normalization, and commutativity failures.
Different timing case counts reflect the growing `2N^2` exact expansion and are
not throughput comparisons.

## Safe Float64x5-x8 addition — 2026-08-12 — Zen 3

| Width | Max observed `|err|/(u^N(|x|+|y|))` |
|---|---:|
| x5 | 0.0528336 |
| x6 | 0.0254805 |
| x7 | 0.0109878 |
| x8 | 0.00541727 |

All measured widths had zero oracle, normalization, and commutativity failures.

## Existing FMA context

The x2-x4 fixed-cost/safe FMA snapshots remain architecture-sensitive: scalar x3
and safe x4 are not current performance wins, while SIMD x3/x4 paths can be
faster. The old x4 FastTwoSum network remains rejected due its concrete
cancellation defect.

## Interpretation

- direct higher-limb FMA is numerically distinct from rounded mul-then-add;
- cancellation demonstrates that distinction most strongly;
- safe x5-x8 add/mul and x5 FMA are correctness baselines, not optimized APIs;
- each width/operation keeps a separate empirical error metric;
- formal tail analysis precedes fixed-cost minimization;
- deployment CPU and downstream solver A/B remain required for performance
  conclusions.
