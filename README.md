# MultiFloatArithmetic.jl

[![CI](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml)

Verification-oriented arithmetic research kernels for
[`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The package explores fixed-length arithmetic for a future MultiFloat-native
linear algebra backend and downstream solvers such as SDPX.jl. Correctness gates
come before instruction-count wins.

## Current stage

M3 safe addition is frozen for Float64x5-x8 and M4 now has its first safe
multiplication baseline for **Float64x5**.

- x2 has a pinned structural FMA proof;
- x3 has a proof-driven fixed-cost normalization repair and structural proof;
- x4 uses a conservative TwoSum + `renormalize` baseline after a concrete
  FastTwoSum cancellation defect was found;
- x5-x8 have exact-rational Experimental add/sub/mul/FMA oracles;
- x5-x8 share one no-FastTwoSum safe addition implementation;
- Float64x5 has an over-complete, exact-tail-accounted, commutative safe
  multiplication baseline;
- formal higher-limb bounds, optimized multiplication/FMA, and native linear
  algebra remain future work.

See [STATUS.md](STATUS.md), [docs/REFERENCE_X5_X8.md](docs/REFERENCE_X5_X8.md),
[docs/ADD_SAFE_X5_X8.md](docs/ADD_SAFE_X5_X8.md), and
[docs/MUL5_SAFE.md](docs/MUL5_SAFE.md).

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

## Safe Float64x5-x8 addition baseline

```julia
safe_sum = E.add_safe(x8, y8)
# or width-specific wrappers E.add5_safe ... E.add8_safe
```

For width N, `add_safe` uses N general TwoSums, fully renormalizes the exact 2N
terms, keeps the N-limb head, and renormalizes it. No FastTwoSum assumption is
introduced.

The permanent corpora require exact head+discarded-tail accounting, normalized
output, bitwise commutativity, and bitwise equality with `reference_add`.
First observed operand-relative constants were approximately 0.05283 (x5),
0.02548 (x6), 0.01099 (x7), and 0.00542 (x8) in units of
`u^N(|x|+|y|)`. CI uses C=1 only as an empirical regression gate.

## Safe Float64x5 multiplication baseline

```julia
T5 = MultiFloat{Float64,5}
x5 = T5(BigFloat("1.2345678901234567890123456789012345"))
y5 = T5(BigFloat("0.8765432109876543210987654321098765"))
product5 = E.mul5_safe(x5, y5)
```

`mul5_safe` evaluates all 25 limb pairs with `two_prod`, keeps both product and
residual, canonicalizes the resulting 50 finite terms so operand swap is
bitwise deterministic, fully renormalizes, truncates to five limbs, and
renormalizes the head. It does not use FastTwoSum.

The permanent corpus checks every individual TwoProd exactly, then requires
exact five-limb-head + discarded-tail reconstruction of `x*y`, normalized
output, bitwise commutativity, and bitwise equality with `reference_mul`.
The first diagnostic measured a worst observed `|err|/(u^5|x*y|)` of about
**0.04841**, with zero oracle/normalization/commutativity failures. CI uses C=1
as an empirical gate, not a theorem.

The first Zen 3 scalar timing was about **33.6 ms for 500 products**. This is
intentionally a correctness baseline, not a performance kernel. Its current
contract also excludes pair-product underflow/subnormal TwoProd components until
a dedicated underflow proof is available.

## Experimental namespace

Experimental APIs include rejected performance candidates, exact correctness
oracles, and safe higher-limb baselines. They may change without deprecation
during the 0.x series and are not production performance recommendations.

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
include("benchmark/division_digits.jl")
include("benchmark/codegen.jl")
```

## Roadmap

1. Derive formal discarded-tail bounds for the safe addition family before
   fixed-cost addition minimization.
2. **M4:** generalize the safe multiplication baseline from x5 through x8,
   measure each width independently, then derive a formal multiplication bound
   before fixed-cost search.
3. **M5:** direct FMA/submul x5 through x8.
4. **M6:** reciprocal/division/sqrt via precision doubling.
5. **M7:** MultiFloatVec-native DOT/SYRK/TRSM/GEMM.
6. Integrate only accepted arithmetic into downstream solvers.

Research software; passing CI is reproducible evidence, not a production/formal
guarantee.

MIT.
