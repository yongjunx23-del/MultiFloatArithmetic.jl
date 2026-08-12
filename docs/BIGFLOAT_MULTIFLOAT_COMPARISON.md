# BigFloat versus MultiFloat precision/speed comparison

This benchmark compares the user-requested pairs:

| Pair | BigFloat nominal precision | MultiFloat nominal significand capacity |
|---|---:|---:|
| BigFloat128 vs Float64x2 | 128 bits | about 106 bits |
| BigFloat256 vs Float64x4 | 256 bits | about 212 bits |
| BigFloat512 vs Float64x8 | 512 bits | about 424 bits |

These are intentionally **not equal-precision pairs**. Each BigFloat setting has
about 20.8% more nominal significand bits than the corresponding Float64xN.
Therefore a speed result must be read together with the precision result.

## Implementations compared

For Float64x2 and Float64x4, `+`, `*`, and `x*y+c` use the original
MultiFloats.jl arithmetic implementation. The benchmark additionally records this
project's `fma_fast` timing and precision for x2/x4.

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
not a formal worst-case bound.

## Speed methodology

Inputs are converted before timing. The benchmark reports the best of several
hosted-runner scalar batches in nanoseconds per operation. BigFloat arithmetic is
run under the requested `setprecision` value. MultiFloat x2/x4 uses upstream
operators directly; x8 uses the project's safe baselines.

Hosted-runner timing is informational only. CPU microarchitecture, Julia/MPFR
versions, allocation behavior, and compiler state can materially change ratios.
The CI log is the authoritative snapshot for each commit.

Reproduce with:

```julia
include("benchmark/bigfloat_multifloat_comparison.jl")
```
