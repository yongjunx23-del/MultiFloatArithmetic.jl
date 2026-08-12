# MultiFloatArithmetic.jl

[![CI](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml)

Verification-oriented arithmetic research kernels for
[`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The package explores fixed-length arithmetic for future MultiFloat-native linear
algebra and downstream solvers such as SDPX.jl. Correctness gates come before
instruction-count wins.

## Current stage

M3/M4/M5 correctness baselines now cover **Float64x5 through Float64x8** for
addition, multiplication, and direct FMA.

- x2/x3 have fixed-cost FMA structural verification;
- x4 uses a cancellation-safe FMA fallback after a concrete FastTwoSum defect;
- x5-x8 have exact-rational Experimental add/sub/mul/FMA oracles;
- `add_safe` and `mul_safe` provide no-FastTwoSum x5-x8 correctness families;
- `fma_safe` forms one exact direct `x*y+c` expansion before N-limb truncation
  for N=5:8;
- formal higher-limb tail bounds, fixed-cost high-limb networks, reciprocal/
  division/sqrt, and native linear algebra remain future work.

The current high-limb safe implementations are correctness references, not hot
kernels. In particular, the requested BigFloat comparison shows that the
intentionally over-complete x8 safe multiplication is far slower than
BigFloat512 and must be replaced before M7 performance integration.

See [STATUS.md](STATUS.md), [docs/REFERENCE_X5_X8.md](docs/REFERENCE_X5_X8.md),
[docs/ADD_SAFE_X5_X8.md](docs/ADD_SAFE_X5_X8.md),
[docs/MUL_SAFE_X5_X8.md](docs/MUL_SAFE_X5_X8.md),
[docs/FMA_SAFE_X5_X8.md](docs/FMA_SAFE_X5_X8.md), and
[docs/BIGFLOAT_MULTIFLOAT_COMPARISON.md](docs/BIGFLOAT_MULTIFLOAT_COMPARISON.md).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/yongjunx23-del/MultiFloatArithmetic.jl")
```

## Accepted research FMA API

```julia
using MultiFloats
using MultiFloatArithmetic

x = Float64x4(BigFloat("1.234567890123456789"))
y = Float64x4(BigFloat("0.987654321098765432"))
c = Float64x4(BigFloat("0.125"))
z = fma_fast(x, y, c)
```

`fma_fast` supports 2-, 3-, and 4-limb values. It is an operand-relative research
kernel, not a correctly rounded FMA.

## Higher-limb exact oracle

```julia
const E = MultiFloatArithmetic.Experimental
T8 = MultiFloat{Float64,8}
x8 = T8(BigFloat("1.2345678901234567890123456789"))
y8 = T8(BigFloat("0.9876543210987654321098765432"))
c8 = T8(BigFloat("-0.1250000000000000000000000001"))

oracle_sum = E.reference_add(x8, y8)
oracle_product = E.reference_mul(x8, y8)
oracle_fma = E.reference_fma(x8, y8, c8)
```

`reference_add/sub/mul/fma` support Float32/Float64 N=5:8 and deliberately use
exact `Rational{BigInt}` arithmetic before one final pack.

## Safe higher-limb add/mul

```julia
safe_sum = E.add_safe(x8, y8)
safe_product = E.mul_safe(x8, y8)
```

Both families preserve and normalize over-complete exact expansions before
N-limb truncation and introduce no FastTwoSum assumptions. They are correctness
baselines, not hot APIs.

## Direct Float64x5-x8 FMA baseline

```julia
safe_fma = E.fma_safe(x8, y8, c8)
# Width-specific wrappers E.fma5_safe ... E.fma8_safe are also available.
```

For width N, `fma_safe` keeps all `2N^2` exact TwoProd product/residual components
and all N limbs of `c`, canonicalizes the resulting `2N^2+N` component direct
`x*y+c` expansion, fully renormalizes it, and only then truncates to N limbs.
It does not call the rounded N-limb `mul_safe` result and does not use FastTwoSum.

The permanent corpus requires exact full-expansion head+tail accounting,
normalized output, bitwise x/y symmetry, and bitwise equality with
`reference_fma`, including deep destructive cancellation.

First Zen 3 maximum observed `|err|/(u^N(|xy|+|c|))` values were approximately:

- x5: 0.04928
- x6: 0.02276
- x7: 0.00452
- x8: 0.000826

The reason for keeping a direct FMA is empirical and strong: the rounded
composition `add_safe(mul_safe(x,y),c)` disagreed with the direct exact-FMA oracle
in every destructive-cancellation case tested at x5, x6, x7, and x8, while the
direct path had zero oracle mismatches.

The safe direct implementation is still slower than safe mul-then-add. M5 freezes
correctness first; fixed-cost fused optimization follows formal tail/error
analysis.

## BigFloat versus MultiFloat comparison

A permanent benchmark compares the requested practical pairs:

- BigFloat128 vs Float64x2 (~106 nominal significand bits),
- BigFloat256 vs Float64x4 (~212 bits),
- BigFloat512 vs Float64x8 (~424 bits).

The BigFloat settings have about 20.8% more nominal precision, so speed and
precision are reported together. On the first Zen 3 snapshot, original upstream
MultiFloats x2/x4 arithmetic had large throughput advantages over BigFloat128/256,
while the current correctness-only x8 safe multiplication was about 3,600x
slower than BigFloat512. See the comparison document for methodology and exact
numbers.

## Experimental namespace

Experimental APIs include rejected performance candidates, exact correctness
oracles, and safe higher-limb baselines. They may change without deprecation and
are not production performance recommendations.

## Reproducing evidence

```julia
using Pkg
Pkg.test("MultiFloatArithmetic")

include("benchmark/smoke.jl")
include("benchmark/simd_widths.jl")
include("benchmark/fma3_repair.jl")
include("benchmark/fma4_safe.jl")
include("benchmark/add5_safe.jl")
include("benchmark/add6_add8_safe.jl")
include("benchmark/mul5_safe.jl")
include("benchmark/mul6_mul8_safe.jl")
include("benchmark/fma5_safe.jl")
include("benchmark/fma6_fma8_safe.jl")
include("benchmark/bigfloat_multifloat_comparison.jl")
include("benchmark/division_digits.jl")
include("benchmark/codegen.jl")
```

## Roadmap

1. Derive formal discarded-tail/error bounds for safe x5-x8 add/mul/FMA and use
   them to drive smaller verified fixed-cost networks, with x8 multiplication a
   priority because the current correctness baseline is extremely slow.
2. **M6:** reciprocal/division/sqrt correctness baselines via precision doubling
   and direct residual correction where useful.
3. **M7:** MultiFloatVec-native DOT/SYRK/TRSM/GEMM and downstream solver A/B only
   after competitive arithmetic survivors exist.

Research software; passing CI is reproducible evidence, not a production/formal
guarantee.

MIT.
