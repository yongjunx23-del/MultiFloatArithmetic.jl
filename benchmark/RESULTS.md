# Recorded benchmark results

These are informational snapshots from GitHub-hosted runners. They are retained
to make architecture dependence visible; they are not hard performance gates.
All ratios use the minimum of seven elapsed-time samples after warm-up.

## Run 21 — 2026-08-11 — Ice Lake

Julia 1.10.11, `Sys.CPU_NAME = icelake-client`.

### FMA: `mul+add/fused`

| Type | Scalar | Vec4 smoke | Vec2 | Vec4 | Vec8 |
|---|---:|---:|---:|---:|---:|
| Float64x2 | 1.087 | 1.077 | 1.071 | 1.085 | 1.033 |
| Float64x3 | 0.831 | 1.132 | 1.255 | 1.063 | 1.080 |
| Float64x4 | 0.900 | 1.260 | 2.017 | 1.157 | 1.155 |

Values above one favor `fma_fast`.

### Quotient-digit division

| Type | Scalar upstream/specialized | Vec4 upstream/specialized | Full/specialized scalar | Full/specialized Vec4 |
|---|---:|---:|---:|---:|
| Float64x2 | 0.380 | 0.847 | 0.981 | 1.006 |
| Float64x3 | 0.429 | 0.918 | 1.026 | 1.005 |
| Float64x4 | 0.237 | 0.683 | 1.042 | 1.008 |

For `upstream/specialized`, values below one mean upstream `/` is faster. For
`full/specialized`, values above one favor the specialized residual product.

## Run 24 — 2026-08-12 — Zen 3

PR #7, Julia 1.10.11, `Sys.CPU_NAME = znver3`.

### FMA: `mul+add/fused`

| Type | Scalar | Vec4 smoke | Vec2 | Vec4 | Vec8 |
|---|---:|---:|---:|---:|---:|
| Float64x2 | 1.020 | 1.005 | 1.252 | 1.051 | 0.942 |
| Float64x3 | 0.488 | 1.891 | 2.197 | 1.972 | 1.592 |
| Float64x4 | 0.510 | 2.529 | 2.537 | 2.537 | 1.950 |

The x2 Vec8 result is a small regression; the other vector cases favor the
candidate. Scalar x3/x4 regress strongly.

### Quotient-digit division

| Type | Scalar upstream/specialized | Vec4 upstream/specialized | Full/specialized scalar | Full/specialized Vec4 |
|---|---:|---:|---:|---:|
| Float64x2 | 0.163 | 0.880 | 1.015 | 1.018 |
| Float64x3 | 0.195 | 1.366 | 1.007 | 1.008 |
| Float64x4 | 0.126 | 0.848 | 1.096 | 1.013 |

The specialized quotient candidate loses badly for every scalar type. Vec4 is
mixed: x2/x4 lose, while x3 wins on this runner. That inconsistency supports
keeping it experimental rather than replacing upstream division.

## Interpretation

The stable conclusions across both snapshots are narrower than either single
run:

- scalar x3/x4 FMA is not a performance win;
- vector FMA is promising, especially x3/x4, but must be selected by target and
  width rather than assumed universally faster;
- quotient-digit scalar division is not competitive;
- quotient-digit Vec4 performance is architecture/type dependent;
- lower instruction count alone does not establish lower wall time.

Before downstream use, rerun on the deployment CPU and perform an end-to-end
solver A/B with iterations, residuals, certificates, time, and memory.
