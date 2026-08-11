# MultiFloatArithmetic.jl

Experimental, verification-oriented arithmetic kernels for [`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The project is intentionally small and research-first. Its immediate goal is to evaluate branch-free fused multiply-add (FMA/MAC) arithmetic for `Float64x2`, `Float64x3`, and `Float64x4`, then use the same methodology to explore provably-correct higher-limb arithmetic and a dedicated MultiFloat linear algebra backend.

## Current scope

- `Float64x2`, `Float64x3`, `Float64x4` fused multiply-add research kernels.
- Scalar and `MultiFloatVec` execution through the same arithmetic network.
- BigFloat differential testing.
- Commutativity and SIMD lane-equivalence checks.
- Explicit destructive-cancellation tests.
- Lightweight hosted-runner timing for `fma_fast(x,y,c)` versus `x*y+c`.

## Numerical contract

The current `fma_fast` path is an **operand-relative** fast FMA research kernel. It must not be interpreted as a correctly-rounded or strong result-relative FMA under severe cancellation. Cancellation-sensitive residual, refinement, and certification code should keep a stronger arithmetic path unless and until a stronger formal contract is established.

The arithmetic network must not be algebraically reordered or silently contracted by the compiler without re-verification.

## Roadmap

1. Validate and benchmark x2/x3/x4 fused arithmetic.
2. Build counterexample-guided and verifier-compatible search tooling for x5-x8 arithmetic.
3. Develop `MultiFloatVec`-native DOT/SYRK/TRSM/GEMM kernels.
4. Build a standalone `MultiFloatLinearAlgebra.jl`-style backend.
5. Integrate only validated kernels into downstream solvers such as SDPX.jl.

## Status

Research software. No arithmetic kernel should be treated as production-ready solely because randomized tests pass; formal error-bound/non-overlap verification and adversarial testing are separate acceptance gates.

## License

MIT.
