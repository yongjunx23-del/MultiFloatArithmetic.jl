# MultiFloatArithmetic.jl

[![CI](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml)

Verification-oriented arithmetic research kernels for
[`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The package explores fixed-length arithmetic for a future MultiFloat-native
linear algebra backend and downstream solvers such as SDPX.jl. Correctness gates
come before instruction-count wins.

## Current stage

M3/M4 correctness baselines now cover **Float64x5 through Float64x8** for
addition and multiplication.

- x2/x3 have fixed-cost FMA structural verification;
- x4 uses a cancellation-safe FMA fallback after a concrete FastTwoSum defect;
- x5-x8 have exact-rational Experimental add/sub/mul/FMA oracles;
- `add_safe` provides one no-FastTwoSum x5-x8 addition family;
- `mul_safe` provides one no-FastTwoSum x5-x8 multiplication family;
- direct x5-x8 FMA, formal high-limb tail bounds, and native linear algebra are
  the next major stages.

See [STATUS.md](STATUS.md), [docs/REFERENCE_X5_X8.md](docs/REFERENCE_X5_X8.md),
[docs/ADD_SAFE_X5_X8.md](docs/ADD_SAFE_X5_X8.md), and
[docs/MUL_SAFE_X5_X8.md](docs/MUL_SAFE_X5_X8.md).

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
kernel, not a correctly rounded FMA; read `docs/NUMERICAL_CONTRACT.md` before
using it in residual/certificate code.

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

## Safe Float64x5-x8 addition

```julia
safe_sum = E.add_safe(x8, y8)
```

For width N, `add_safe` uses N general TwoSums, fully renormalizes the exact 2N
terms, keeps the N-limb head, and renormalizes it. No FastTwoSum assumption is
introduced. First observed operand-relative constants were approximately
0.05283, 0.02548, 0.01099, and 0.00542 for x5-x8.

## Safe Float64x5-x8 multiplication

```julia
safe_product = E.mul_safe(x8, y8)
# or E.mul5_safe ... E.mul8_safe
```

For width N, `mul_safe` evaluates every N×N limb pair with `two_prod`, preserves
all product/residual components (`2N^2` terms), canonicalizes their order under
operand swap, fully renormalizes, keeps the N-limb head, and renormalizes it.
No FastTwoSum is used.

The permanent corpus checks every TwoProd exactly, exact head+discarded-tail
reconstruction of `x*y`, normalized output, bitwise commutativity, and bitwise
`reference_mul` equality. First maximum observed `|err|/(u^N|x*y|)` values were:

- x5: 0.0484069
- x6: 0.0183703
- x7: 0.00725543
- x8: 0.00253803

CI uses C=1 per width as an empirical regression gate, not a theorem. The safe
multiplication family is deliberately slow and currently excludes subnormal or
underflowing TwoProd components pending a dedicated proof.

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
include("benchmark/division_digits.jl")
include("benchmark/codegen.jl")
```

## Roadmap

1. Derive formal discarded-tail bounds for safe x5-x8 add/mul before fixed-cost
   minimization.
2. **M5:** build direct correctness-first FMA/submul x5 through x8, then optimize
   accepted survivors.
3. **M6:** reciprocal/division/sqrt via precision doubling.
4. **M7:** MultiFloatVec-native DOT/SYRK/TRSM/GEMM and downstream solver A/B.

Research software; passing CI is reproducible evidence, not a production/formal
guarantee.

MIT.
