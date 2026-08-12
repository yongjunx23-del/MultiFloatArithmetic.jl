# Recorded benchmark results

These are informational GitHub-hosted-runner snapshots, not portable hard gates.
Ratios use the minimum observed elapsed time after warm-up.

## Earlier FMA snapshots

### 2026-08-11 — Ice Lake

| Type | Scalar | Vec2 | Vec4 | Vec8 |
|---|---:|---:|---:|---:|
| Float64x2 | 1.087x | 1.071x | 1.085x | 1.033x |
| Float64x3 | 0.831x | 1.255x | 1.063x | 1.080x |
| Float64x4 old candidate | 0.900x | 2.017x | 1.157x | 1.155x |

### 2026-08-12 — Zen 3

| Type | Scalar | Vec2 | Vec4 | Vec8 |
|---|---:|---:|---:|---:|
| Float64x2 | 1.020x | 1.252x | 1.051x | 0.942x |
| Float64x3 | 0.488x | 2.197x | 1.972x | 1.592x |
| Float64x4 old candidate | 0.510x | 2.537x | 2.537x | 1.950x |

Values above one favor `fma_fast` over upstream `x*y+c`. These old x4 numbers
are performance evidence only; the old x4 end network is no longer acceptable
because its cancellation FastTwoSum precondition is invalid.

## 2026-08-12 — safe x4 baseline — generic hosted CPU

Julia 1.10.11.

### Safe FMA versus upstream

| Type | Scalar | Vec2 | Vec4 | Vec8 |
|---|---:|---:|---:|---:|
| Float64x2 | 1.008x | 1.135x | 1.030x | 1.017x |
| Float64x3 repaired | 0.333x | 1.768x | 1.417x | 1.067x |
| Float64x4 safe | 0.284x | 1.591x | 1.586x | 1.311x |

### x4 old versus safe baseline

| Workload | Old | Safe | safe/old |
|---|---:|---:|---:|
| scalar, 20k | 0.062 ms | 0.345 ms | 5.593x |
| Vec2, 20k lanes | 0.155 ms | 0.211 ms | 1.364x |
| Vec4, 20k lanes | 0.080 ms | 0.109 ms | 1.352x |
| Vec8, 20k lanes | 0.064 ms | 0.068 ms | 1.056x |

Ordinary random inputs produced no bitwise old/safe differences in 20,000
scalar samples or the SIMD samples.

### Exact-cancellation diagnostic

For 10,000 seeded Float64x4 inputs with `c = Float64x4(-big(x)*big(y))`:

- old first-FastTwoSum precondition violated: **3,755**;
- old versus safe outputs changed: **5,258**;
- old non-normalized outputs: **5,258**;
- safe non-normalized outputs: **0**.

This corpus is the strongest current evidence for keeping the correctness-first
x4 fallback despite its scalar cost.

## Quotient-digit division

Across recorded architectures, the specialized quotient-digit scalar path is
consistently much slower than upstream division. Vec4 results are mixed by type
and architecture. It remains `Experimental`.

## Interpretation

Stable conclusions:

- scalar x3/x4 custom FMA is not a performance win;
- x2 scalar benefit is marginal and architecture-sensitive;
- SIMD x3 and safe x4 remain promising, but selection must be target/width
  specific;
- instruction count alone is not sufficient evidence;
- cancellation correctness overrides the larger speedups measured for the old
  x4 candidate.

Before downstream use, rerun on the deployment CPU and perform an end-to-end
solver A/B including iterations, residuals, certificates, time, and memory.
