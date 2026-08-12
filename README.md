# MultiFloatArithmetic.jl

[![CI](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml)

Verification-oriented arithmetic research kernels for
[`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The package is focused on fixed-length, branch-free arithmetic that can later
support a MultiFloat-native linear algebra backend and downstream solvers such
as SDPX.jl. It is deliberately conservative: a kernel is not promoted merely
because randomized tests pass or its instruction count is smaller.

## Current stage

The project is at **M1: empirical validation and API freeze for x2/x3/x4 fused
multiply-add**. The x5-x8 reference path, formal verification, higher-limb
multiplication/division, and native linear algebra backend have not yet been
implemented.

Two hosted-runner studies on different CPU targets produced these decisions:

- keep the x2 scalar FMA as a marginal, architecture-dependent candidate;
- keep x3/x4 scalar FMA opt-in because they regressed on both runners;
- continue SIMD FMA as the main optimization path, but do not auto-select it:
  the tested vector cases were mostly faster, with one slight x2 Vec8 regression;
- retain quotient-digit division only under `Experimental`: its scalar path was
  consistently much slower, while Vec4 results were architecture-dependent.

See [STATUS.md](STATUS.md) for the acceptance state and
[benchmark/RESULTS.md](benchmark/RESULTS.md) for the recorded measurements.

## Installation

The package is not registered. Install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/yongjunx23-del/MultiFloatArithmetic.jl")
```

## Accepted top-level API

```julia
using MultiFloats
using MultiFloatArithmetic

x = Float64x4(BigFloat("1.234567890123456789"))
y = Float64x4(BigFloat("0.987654321098765432"))
c = Float64x4(BigFloat("0.125"))

z = fma_fast(x, y, c)
```

`fma_fast` supports `Float64x2`, `Float64x3`, `Float64x4`, and matching
`MultiFloatVec` values. The implementation is generic enough for binary32 and
has smoke coverage there, but the present optimization target and benchmark
decisions are binary64.

## Numerical contract

`fma_fast` is an **operand-relative** FMA research kernel. It is not a
correctly-rounded FMA and does not provide a strong result-relative guarantee
under destructive cancellation. Do not use it by default in residual,
refinement, stopping-criterion, or certification code.

The arithmetic network must not be algebraically reordered, contracted beyond
the explicit error-free transforms, or compiled under `@fastmath` without
re-verification. Read [docs/NUMERICAL_CONTRACT.md](docs/NUMERICAL_CONTRACT.md)
before downstream integration.

## Experimental namespace

Rejected or not-yet-accepted candidates remain available only for reproducible
research:

```julia
using MultiFloats
using MultiFloatArithmetic

x = Float64x4(1.25)
y = Float64x4(0.75)

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
include("benchmark/division_digits.jl")
include("benchmark/codegen.jl")
```

Hosted-runner timings are informational. A downstream kernel should be accepted
only after architecture-specific repeated measurements and an end-to-end solver
A/B.

## Roadmap

1. Freeze x2/x3/x4 behavior and preserve a permanent adversarial corpus.
2. Translate the accumulation networks into verifier inputs and obtain formal
   non-overlap/error-bound evidence.
3. Build a correctness-first x5 reference path before optimizing x5-x8.
4. Add verified multiplication and direct FMA incrementally from x5 to x8.
5. Implement reciprocal/division/square root through precision doubling.
6. Build `MultiFloatVec`-native DOT/SYRK/TRSM/GEMM kernels.
7. Integrate only accepted kernels into downstream solvers.

## Status and license

This remains research software. Passing CI establishes reproducible empirical
evidence, not a formal proof or production guarantee.

MIT.
