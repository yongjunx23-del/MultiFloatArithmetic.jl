# MultiFloatArithmetic.jl

[![CI](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml)

Verification-oriented arithmetic research kernels for
[`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The package explores fixed-length arithmetic for a future MultiFloat-native
linear algebra backend and downstream solvers such as SDPX.jl. Correctness gates
come before instruction-count wins.

## Current stage

The project is at **M1.1: structural verification and cancellation hardening for
x2/x3/x4 FMA**.

- x2 has a pinned structural proof;
- x3 has a proof-driven fixed-cost normalization repair and pinned structural
  proof;
- the original x4 fixed-cost end network had an invalid FastTwoSum assumption
  under cancellation;
- x4 currently uses a conservative TwoSum + `renormalize` baseline while a
  source-specific fixed-cost replacement remains research;
- x5-x8, formal error constants, and native linear algebra are not implemented.

See [STATUS.md](STATUS.md) for the acceptance state and
[benchmark/RESULTS.md](benchmark/RESULTS.md) for measurements.

## Installation

The package is not registered:

```julia
using Pkg
Pkg.add(url="https://github.com/yongjunx23-del/MultiFloatArithmetic.jl")
```

## API

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

Rejected or not-yet-accepted candidates remain available only for reproducible
research:

```julia
z = MultiFloatArithmetic.Experimental.div_digits(x, y)
```

The `Experimental` API may change without deprecation during the 0.x series and
must not be treated as the package's performance recommendation.

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

1. Freeze x2/x3 and the cancellation-safe x4 baseline with permanent adversarial
   corpora.
2. Prove the current empirical error constants or replace them with proved ones.
3. Derive a source-specific fixed-cost x4 compression if it beats the safe
   baseline without reintroducing hidden magnitude assumptions.
4. Build a correctness-first x5 reference implementation, then extend toward
   x6-x8.
5. Add reciprocal/division/square-root strategies at higher limb counts.
6. Build MultiFloatVec-native DOT/SYRK/TRSM/GEMM kernels.
7. Integrate only accepted arithmetic into downstream solvers.

## Status and license

Research software; passing CI is reproducible empirical evidence, not a formal
production guarantee.

MIT.
