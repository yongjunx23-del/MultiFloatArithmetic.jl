# Recorded benchmark results

These are informational GitHub-hosted-runner snapshots, not portable hard gates.

## M7.2 blocked Float64x4 Cholesky — 2026-08-12 — Zen 3

The optimized lower-factorization schedule preserves the scalar baseline's exact
per-entry FMA chronology and operand roles. Before timing, every tested result
must be whole-matrix bitwise `===` to the unblocked kernel and normalized.

The public route uses Vec8, block size 8, only for dense `Matrix{Float64x4}` lower
factorization with `n >= 32`; all other cases fall back to the unblocked kernel.

| n | public POTRF | unblocked | speedup |
|---:|---:|---:|---:|
| 32 | 0.212 ms | 0.365 ms | **1.72x** |
| 64 | 1.058 ms | 2.468 ms | **2.33x** |
| 96 | 2.850 ms | 7.801 ms | **2.74x** |

A previous candidate run tested BS=8/16/24/32; BS=8 was consistently fastest.
Permanent n=17/33 regressions, an n=73 remainder case and all benchmarked block
sizes remained bitwise equal to unblocked. Linux Julia 1.10/current and current
macOS passed with the public optimized dispatch enabled.

## M7.1 dense GEMM microkernel — 2026-08-12 — Zen 3

The public Float64x4 dense GEMM route uses an explicit `MultiFloatVec` MR=8,
NR=2 microkernel for unit alpha and sufficiently large dense matrices. Its
per-output k reduction remains ordered, and every candidate output is required
to be bitwise equal to the streaming baseline.

| n | public x4 GEMM | streaming | public speedup | generic `mul!` | public/generic |
|---:|---:|---:|---:|---:|---:|
| 16 | 0.038 ms | 0.045 ms | 1.19x | 0.290 ms | 7.66x |
| 32 | 0.277 ms | 0.347 ms | 1.25x | 2.351 ms | 8.48x |
| 48 | 0.923 ms | 1.188 ms | 1.29x | 7.890 ms | 8.55x |
| 64 | 2.205 ms | 2.796 ms | 1.27x | 18.653 ms | 8.46x |

Explicit Vec4 was rejected. Explicit x2 vector packing was not promoted because
the compiler-generated streaming path was already strong; in the same run x2
public GEMM was about 7.3-9.3x faster than generic `mul!` across n=16:64.

Other M7 snapshot results from the same runner:

- Float64x2 DOT: 0.153 ms vs 0.253 ms generic -> 1.65x;
- Float64x4 DOT: 0.927 ms vs 1.577 ms generic -> 1.70x;
- Float64x2 GEMV 256x128: 0.039 ms vs 0.066 ms -> 1.71x;
- Float64x4 GEMV 256x128: 0.344 ms vs 0.608 ms -> 1.77x.

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

The current five-pass width-specialized x4 FMA supersedes the earlier
pre-audit x4 timing recorded in the original BigFloat comparison. On a current
Zen 3 run, project x4 direct FMA was ~10.8 ns/op versus ~18.8 ns/op for upstream
`x*y+c`; the requested BigFloat256 mul+add was ~155.5 ns/op in that run. Keep
historical snapshots tied to their exact code revision rather than mixing them.

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

## Interpretation

- direct higher-limb FMA is numerically distinct from rounded mul-then-add;
- cancellation demonstrates that distinction at every x5-x8 width;
- x4 five-pass direct FMA is now a competitive scalar/SIMD primitive;
- matrix-level ordered FMA scheduling produces much larger wins than scalar
  add/mul measurements alone suggest;
- explicit SIMD is accepted only after bitwise-equality A/B and stable benefit;
- current safe x8 add/mul/FMA is correctness-only and is not a competitive hot
  path, especially multiplication;
- deployment CPU and downstream solver A/B remain required before performance
  claims are generalized beyond the recorded hosts.
