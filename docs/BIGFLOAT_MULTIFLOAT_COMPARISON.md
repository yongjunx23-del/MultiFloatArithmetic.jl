# BigFloat versus MultiFloat precision/speed comparison

This benchmark compares the user-requested pairs:

| Pair | BigFloat nominal precision | MultiFloat nominal significand capacity |
|---|---:|---:|
| BigFloat128 vs Float64x2 | 128 bits | about 106 bits |
| BigFloat256 vs Float64x4 | 256 bits | about 212 bits |
| BigFloat512 vs Float64x8 | 512 bits | about 424 bits |

These are intentionally **not equal-precision pairs**. Each BigFloat setting has
about 20.8% more nominal significand bits than the corresponding Float64xN.
Therefore speed and precision must be read together.

## Implementations compared

For Float64x2 and Float64x4, `+`, `*`, and `x*y+c` use the original
MultiFloats.jl v3.2.6 arithmetic implementation. The benchmark additionally
records this project's `fma_fast` timing and precision for x2/x4.

Upstream MultiFloats.jl does not provide arithmetic for Float64x8. The x8 rows
therefore use this project's correctness-first `Experimental.add8_safe` and
`Experimental.mul8_safe` baselines. They are intentionally branchy/allocating
and are not performance kernels.

## Precision methodology

A fixed seeded corpus of exact positive `Rational{BigInt}` source values is
shared by both formats. Each source is independently rounded into the requested
BigFloat precision and into Float64xN. Addition, multiplication, and (where
available) multiply-plus-add are then compared end-to-end against the exact
rational source expression.

The benchmark reports maximum observed relative error and
`-log2(max relative error)` as an observed effective-bit indicator. This metric
includes both input representation error and arithmetic error. It is empirical,
not a formal worst-case bound. Observed effective bits can exceed the simple
`53N` capacity heuristic on a finite corpus and must not be read as a guaranteed
precision theorem.

## First measured snapshot

GitHub hosted runner, Zen 3 (`znver3`), Julia 1.10.11, MultiFloats v3.2.6,
2026-08-12.

### Precision

| Pair / operation | BigFloat max relative error | BigFloat observed bits | MultiFloat max relative error | MultiFloat observed bits |
|---|---:|---:|---:|---:|
| 128 vs x2 add | 3.7989e-39 | 127.6 | 1.5854e-32 | 105.6 |
| 128 vs x2 mul | 5.0334e-39 | 127.2 | 1.6867e-32 | 105.5 |
| 128 vs x2 mul+add | 6.3065e-39 | 126.9 | 1.9117e-32 | 105.4 |
| 256 vs x4 add | 1.3569e-77 | 255.3 | 1.7603e-65 | 215.1 |
| 256 vs x4 mul | 1.7718e-77 | 255.0 | 2.0477e-64 | 211.6 |
| 256 vs x4 mul+add | 1.9863e-77 | 254.8 | 2.0457e-64 | 211.6 |
| 512 vs x8 add | 7.0562e-155 | 512.1 | 8.0722e-131 | 432.2 |
| 512 vs x8 mul | 7.1740e-155 | 512.1 | 1.0228e-130 | 431.8 |

For x2, this project's `fma_fast` measured `2.1053e-32` maximum relative error
(~105.2 observed bits), versus `1.9117e-32` (~105.4 bits) for original
MultiFloats `x*y+c`.

For x4, `fma_fast` measured `2.5128e-64` (~211.3 bits), versus
`2.0457e-64` (~211.6 bits) for original MultiFloats `x*y+c`.

The BigFloat rows are more accurate here largely because the requested BigFloat
precisions contain about 20.8% more nominal significand bits. This comparison is
therefore a requested practical pairing, not an equal-bit algorithm contest.

### Speed

The script times batched independent scalar operations and reports the best
**throughput-equivalent ns/op**. These numbers are not single-operation latency;
sub-nanosecond x2 values can arise from compiler scheduling/vectorization across
independent loop iterations.

| Pair / operation | BigFloat | MultiFloat | Interpretation |
|---|---:|---:|---|
| 128 vs x2 add | 25.8 ns/op | 0.7 ns/op original | original x2 ~35.9x higher throughput |
| 128 vs x2 mul | 35.7 ns/op | 0.4 ns/op original | original x2 ~91.4x higher throughput |
| 128 vs x2 mul+add | 61.4 ns/op | 1.0 ns/op original | original x2 ~61x higher throughput |
| 256 vs x4 add | 33.5 ns/op | 5.6 ns/op original | original x4 ~6.0x higher throughput |
| 256 vs x4 mul | 49.8 ns/op | 5.6 ns/op original | original x4 ~8.9x higher throughput |
| 256 vs x4 mul+add | 80.5 ns/op | 10.7 ns/op original | original x4 ~7.5x higher throughput |
| 512 vs x8 add | 42.5 ns/op | 194.8 ns/op safe baseline | BigFloat512 ~4.6x higher throughput |
| 512 vs x8 mul | 81.4 ns/op | 296637.9 ns/op safe baseline | BigFloat512 ~3644x higher throughput |

Project `fma_fast` throughput-equivalent timing was 0.6 ns/op at x2 and
28.3 ns/op at x4. Thus it beats original x2 `x*y+c` in this batched scalar
snapshot, but is about 2.6x slower than original x4 `x*y+c`, consistent with the
project's earlier conclusion that scalar x4 is not currently a speed path.

The x8 result is the most actionable engineering signal. `mul8_safe` constructs,
sorts, and repeatedly renormalizes all 128 exact TwoProd components. Its current
~296.6 microseconds/product is a correctness reference implementation, not a
competitive arithmetic kernel. **No x8 performance claim should be made until a
smaller verified multiplication network replaces this baseline.**

## Speed methodology

Inputs are converted before timing. BigFloat arithmetic runs under the requested
`setprecision` value. MultiFloat x2/x4 uses upstream operators directly; x8 uses
the project's safe baselines. The benchmark takes the best of several batches,
so the result is a best-case hosted-runner throughput snapshot.

CPU microarchitecture, Julia/MPFR versions, allocation behavior, compiler state,
and SIMD/autovectorization can materially change ratios. Downstream solver
benchmarks remain the final performance criterion.

## Reproduce

```julia
include("benchmark/bigfloat_multifloat_comparison.jl")
```
