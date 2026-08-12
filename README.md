# MultiFloatArithmetic.jl

[![CI](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml)

Verification-oriented arithmetic research kernels for
[`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The package explores fixed-length arithmetic for a future MultiFloat-native
linear algebra backend and downstream solvers such as SDPX.jl. Correctness gates
come before instruction-count wins.

## Current stage

M1.1 is complete and the project is now at **M2: exact x5-x8 reference
arithmetic**.

- x2 has a pinned structural proof;
- x3 has a proof-driven fixed-cost normalization repair and pinned structural
  proof;
- x4 uses a conservative TwoSum + `renormalize` baseline after a concrete
  FastTwoSum cancellation defect was found in the old network;
- x5-x8 now have exact-rational Experimental add/sub/mul/FMA reference oracles;
- optimized x5-x8 networks, formal error constants, and native linear algebra
  are not implemented yet.

See [STATUS.md](STATUS.md) for the acceptance state,
[docs/REFERENCE_X5_X8.md](docs/REFERENCE_X5_X8.md) for the higher-limb oracle,
and [benchmark/RESULTS.md](benchmark/RESULTS.md) for measurements.

## Installation

The package is not registered:

```julia
using Pkg
Pkg.add(url="https://github.com/yongjunx23-del/MultiFloatArithmetic.jl")
```

## Accepted research API

```julia
using MultiFloats
using MultiFloatArithmetic

x = Float64x4(BigFloat("1.234567890123456789"))
y = Float64x4(BigFloat("0.987654321098765432"))
c = Float64x4(BigFloat("0.125"))

z = fma_fast(x, y, c)
```

`fma_fast` supports 2-, 3-, and 4-limb `MultiFloat` values and matching
`MultiFloatVec` values. Binary32 has smoke coverage; binary64 is the primary
optimization target.

## x5-x8 correctness oracle

Higher-limb work begins under `Experimental`, not the accepted hot API:

```julia
using MultiFloats
using MultiFloatArithmetic

const E = MultiFloatArithmetic.Experimental
T = MultiFloat{Float64,8}

x = T(BigFloat("1.2345678901234567890123456789"))
y = T(BigFloat("0.9876543210987654321098765432"))
c = T(BigFloat("0.125"))

s = E.reference_add(x, y)
p = E.reference_mul(x, y)
f = E.reference_fma(x, y, c)
```

`reference_add`, `reference_sub`, `reference_mul`, and `reference_fma` support
Float32/Float64 with 5-8 limbs. They use exact `Rational{BigInt}` arithmetic and
then pack the exact result once. They are intentionally slow and exist only to
validate future optimized kernels. CI cross-checks their packing against an
independent 8192-bit BigFloat result.

## Numerical contract

`fma_fast` is an **operand-relative research FMA**, not a correctly rounded FMA.
It does not promise a strong result-relative guarantee under destructive
cancellation. Read [docs/NUMERICAL_CONTRACT.md](docs/NUMERICAL_CONTRACT.md)
before using it in residual, refinement, stopping, feasibility, or certificate
code.

The x2/x3 arithmetic networks are fixed-cost. The x4 correctness baseline is
intentionally not branch-free because it ends in upstream `renormalize` after a
concrete cancellation defect was found in the old FastTwoSum-based compression.

## Experimental namespace

Experimental APIs include rejected performance candidates and correctness-only
oracles. They may change without deprecation during the 0.x series and must not
be treated as production performance recommendations.

```julia
z = MultiFloatArithmetic.Experimental.div_digits(Float64x4(1.25), Float64x4(0.75))
```

Quotient-digit division remains only for reproducibility; the x5-x8
`reference_*` functions are deliberately slow test/research oracles.

## Reproducing the evidence

```julia
using Pkg
Pkg.test("MultiFloatArithmetic")

include("benchmark/smoke.jl")
include("benchmark/simd_widths.jl")
include("benchmark/fma3_repair.jl")
include("benchmark/fma4_safe.jl")
include("benchmark/division_digits.jl")
include("benchmark/codegen.jl")
```

Hosted-runner timings are informational. A kernel should be selected only after
measurements on the deployment CPU and an end-to-end downstream solver A/B.

## Roadmap

1. **M2:** freeze the exact x5-x8 reference oracle and permanent adversarial
   corpus.
2. **M3:** build safe x5 addition first, measure its discarded tail against the
   oracle, then extend verified addition through x8.
3. **M4:** build commutative multiplication x5 through x8.
4. **M5:** build direct FMA/submul x5 through x8.
5. **M6:** add reciprocal/division/square-root strategies via precision doubling.
6. **M7:** build MultiFloatVec-native DOT/SYRK/TRSM/GEMM kernels.
7. Integrate only accepted arithmetic into downstream solvers.

In parallel, the current x2/x3 empirical error constants still need formal proof
(or replacement with proved constants), and a fixed-cost x4 compression remains
optional research if it can beat the safe baseline without hidden magnitude
assumptions.

## Status and license

Research software; passing CI is reproducible empirical evidence, not a formal
production guarantee.

MIT.
