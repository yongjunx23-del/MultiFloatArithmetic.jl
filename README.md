# MultiFloatArithmetic.jl

[![CI](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml)

Verification-oriented arithmetic research kernels for
[`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The package explores fixed-length arithmetic for future MultiFloat-native linear
algebra and downstream solvers such as SDPX.jl. Correctness gates come before
instruction-count wins.

## Current stage

M3/M4 safe addition and multiplication cover Float64x5-x8. M5 now has the first
**direct Float64x5 FMA correctness baseline**.

- x2/x3 have fixed-cost FMA structural verification;
- x4 uses a cancellation-safe FMA fallback after a concrete FastTwoSum defect;
- x5-x8 have exact-rational Experimental add/sub/mul/FMA oracles;
- `add_safe` and `mul_safe` provide no-FastTwoSum x5-x8 correctness families;
- `fma5_safe` forms one exact 55-component `x*y+c` expansion before five-limb
  truncation;
- safe direct FMA x6-x8, formal high-limb bounds, and native linear algebra are
  the next major stages.

See [STATUS.md](STATUS.md), [docs/REFERENCE_X5_X8.md](docs/REFERENCE_X5_X8.md),
[docs/ADD_SAFE_X5_X8.md](docs/ADD_SAFE_X5_X8.md),
[docs/MUL_SAFE_X5_X8.md](docs/MUL_SAFE_X5_X8.md), and
[docs/FMA5_SAFE.md](docs/FMA5_SAFE.md).

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

oracle_sum = E.reference_add(x8, y8)
oracle_product = E.reference_mul(x8, y8)
```

`reference_add/sub/mul/fma` support Float32/Float64 N=5:8 and deliberately use
exact `Rational{BigInt}` arithmetic before one final pack.

## Safe higher-limb add/mul

```julia
safe_sum = E.add_safe(x8, y8)
safe_product = E.mul_safe(x8, y8)
```

Both families fully preserve/normalize over-complete exact expansions before
N-limb truncation and introduce no FastTwoSum assumptions. They are correctness
baselines, not hot APIs.

## Direct Float64x5 FMA baseline

```julia
T5 = MultiFloat{Float64,5}
x5 = T5(BigFloat("1.234567890123456789012345678901"))
y5 = T5(BigFloat("0.876543210987654321098765432109"))
c5 = T5(BigFloat("-0.5"))

z5 = E.fma5_safe(x5, y5, c5)
```

`fma5_safe` keeps all 50 TwoProd product/residual components from x5
multiplication and all five limbs of `c`, canonicalizes the resulting 55-term
`x*y+c` expansion, fully renormalizes it, and only then truncates to five limbs.
It does not call the rounded five-limb `mul_safe` result and does not use
FastTwoSum.

The permanent corpus requires exact 55-term head+tail accounting, normalized
output, bitwise x/y symmetry, and bitwise equality with `reference_fma`, including
deep destructive cancellation.

First Zen 3 diagnostics showed zero direct-oracle/normalization/symmetry failures
and max observed `|err|/(u^5(|xy|+|c|))` about **0.04928**. By contrast,
`add_safe(mul_safe(x,y),c)` disagreed with the exact direct-FMA oracle in 45/200
ordinary cases, 41/200 scaled cases, and **150/150 destructive-cancellation
cases**. This is why the higher-limb roadmap keeps a direct FMA rather than
assuming separately correct mul/add are an equivalent fused path.

The initial correctness implementation is not faster: 120 direct cases measured
7.707 ms versus 6.231 ms for safe mul-then-add. Fixed-cost direct-FMA optimization
comes after the safe x5-x8 family and formal tail analysis.

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
include("benchmark/division_digits.jl")
include("benchmark/codegen.jl")
```

## Roadmap

1. Derive formal discarded-tail bounds for safe x5-x8 add/mul/FMA before fixed-
   cost minimization.
2. **M5:** generalize direct safe FMA from x5 through x8; then search optimized
   fused survivors.
3. **M6:** reciprocal/division/sqrt via precision doubling.
4. **M7:** MultiFloatVec-native DOT/SYRK/TRSM/GEMM and downstream solver A/B.

Research software; passing CI is reproducible evidence, not a production/formal
guarantee.

MIT.
