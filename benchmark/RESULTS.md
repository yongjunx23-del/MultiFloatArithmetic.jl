# Recorded benchmark results

These are informational GitHub-hosted-runner snapshots, not portable hard gates.

## BigFloat versus MultiFloat — 2026-08-12 — Zen 3

Julia 1.10.11, MultiFloats v3.2.6. The requested pairs are not equal precision:
BigFloat128/256/512 have about 20.8% more nominal significand bits than
Float64x2/x4/x8 (~106/~212/~424 bits).

### End-to-end precision

Maximum relative error against shared exact `Rational{BigInt}` source
expressions:

| Pair / operation | BigFloat observed bits | MultiFloat observed bits |
|---|---:|---:|
| 128 vs x2 add | 127.6 | 105.6 |
| 128 vs x2 mul | 127.2 | 105.5 |
| 128 vs x2 mul+add | 126.9 | 105.4 |
| 256 vs x4 add | 255.3 | 215.1 |
| 256 vs x4 mul | 255.0 | 211.6 |
| 256 vs x4 mul+add | 254.8 | 211.6 |
| 512 vs x8 add | 512.1 | 432.2 |
| 512 vs x8 mul | 512.1 | 431.8 |

The x2/x4 MultiFloat rows use original upstream MultiFloats.jl arithmetic.
At x2, project `fma_fast` observed ~105.2 bits versus ~105.4 for upstream
`x*y+c`; at x4, project `fma_fast` observed ~211.3 bits versus ~211.6 upstream.
These are finite-corpus measurements, not guaranteed precision bounds.

### Batched throughput-equivalent timing

| Pair / operation | BigFloat ns/op | MultiFloat ns/op | Faster side / ratio |
|---|---:|---:|---|
| 128 vs x2 add | 25.8 | 0.7 original | MultiFloat ~35.9x |
| 128 vs x2 mul | 35.7 | 0.4 original | MultiFloat ~91.4x |
| 128 vs x2 mul+add | 61.4 | 1.0 original | MultiFloat ~61x |
| 256 vs x4 add | 33.5 | 5.6 original | MultiFloat ~6.0x |
| 256 vs x4 mul | 49.8 | 5.6 original | MultiFloat ~8.9x |
| 256 vs x4 mul+add | 80.5 | 10.7 original | MultiFloat ~7.5x |
| 512 vs x8 add | 42.5 | 194.8 safe | BigFloat ~4.6x |
| 512 vs x8 mul | 81.4 | 296637.9 safe | BigFloat ~3644x |

These are best batched throughput-equivalent numbers, not single-operation
latencies. In particular, sub-nanosecond x2 values can reflect scheduling or
autovectorization across independent iterations.

Project `fma_fast` measured 0.6 ns/op at x2 and 28.3 ns/op at x4. Thus the x2
kernel is faster than original upstream `x*y+c` in this snapshot, while scalar
x4 is ~2.6x slower than upstream, consistent with prior project conclusions.

Upstream MultiFloats.jl has no x8 arithmetic. The x8 comparison uses this
project's intentionally over-complete safe baselines. The x8 multiplication
result is therefore a **development signal, not a product-speed result**:
`mul8_safe` sorts and renormalizes 128 exact product/residual components and is
orders of magnitude slower than BigFloat512. A smaller verified x8 product
network is required before MultiFloat-native linear algebra can make a credible
x8 performance claim.

See `docs/BIGFLOAT_MULTIFLOAT_COMPARISON.md` for full methodology and raw error
figures.

## Safe Float64x5-x8 direct FMA — 2026-08-12 — Zen 3

The M5 family forms one direct exact `x*y+c` expansion before N-limb truncation.
It had zero direct oracle, normalization, and x/y-symmetry failures at every
width in the first accepted run.

### Direct-FMA error and composition mismatch

| Width | Max observed `|err|/(u^N(|xy|+|c|))` | Ordinary mul-then-add mismatches | Scaled mismatches | Destructive-cancellation mismatches |
|---|---:|---:|---:|---:|
| x5 | 0.0492768 | 45 / 200 | 41 / 200 | **150 / 150** |
| x6 | 0.0227603 | 8 / 50 | 14 / 50 | **48 / 48** |
| x7 | 0.00452122 | 7 / 35 | 9 / 35 | **32 / 32** |
| x8 | 0.000825581 | 6 / 25 | 6 / 25 | **24 / 24** |

`add_safe(mul_safe(x,y),c)` disagreed with the independent exact direct-FMA
oracle in every tested destructive-cancellation case at every width. The direct
path had zero oracle mismatches. This is the central M5 reason to preserve a
true direct FMA for later residual/refinement/certificate work.

CI uses C=1 for the operand-relative metric independently per width as an
empirical regression gate, not a theorem.

### First safe-path timing A/B

| Width | Direct safe FMA | Safe mul-then-add | Composition/direct |
|---|---:|---:|---:|
| x5 | 9.852 ms / 120 | 7.993 ms / 120 | 0.811x |
| x6 | 3.044 ms / 20 | 2.550 ms / 20 | 0.837x |
| x7 | 3.186 ms / 12 | 2.726 ms / 12 | 0.856x |
| x8 | 2.625 ms / 6 | 2.314 ms / 6 | 0.882x |

Different case counts reflect the growing 55/78/105/136-component exact direct
expansions and are not cross-width throughput comparisons. The direct path is
currently slower; it is a correctness reference, not a hot kernel.

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
- cancellation demonstrates that distinction at every x5-x8 width;
- original MultiFloat x2/x4 has a large throughput advantage over the requested
  higher-precision BigFloat pairings in this hosted scalar batch;
- current safe x8 add/mul/FMA is correctness-only and is not a competitive hot
  path, especially multiplication;
- safe x5-x8 add/mul/FMA are oracle-matching baselines for future network
  minimization, not optimized public APIs;
- each width/operation keeps a separate empirical error metric;
- formal tail analysis precedes fixed-cost minimization;
- deployment CPU and downstream solver A/B remain required for performance
  conclusions.
