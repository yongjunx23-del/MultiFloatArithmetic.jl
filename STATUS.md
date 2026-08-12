# Project status

Status date: **2026-08-12**  
Current milestone: **M7 x2-x4 dense linear-algebra baseline; higher-width arithmetic redesign continues**

## Executive assessment

The project has moved beyond an arithmetic-only research package.

Current accepted direction:

- x2/x3 fixed-cost direct-FMA networks with pinned structural verification;
- x4 current five-pass QW direct FMA, reproduced from arXiv:2607.11391 and
  width-specialized for Julia code generation;
- x5-x8 exact-rational oracles and correctness-first add/mul/direct-FMA
  baselines remain Experimental;
- `MFLinearAlgebra` now provides FMA-native DOT/AXPY/GEMV/GEMM/SYRK/TRSV/TRSM
  and Cholesky for Float32/Float64 x2-x4;
- linear-algebra outputs are checked against 1024-bit BigFloat and must remain
  `MultiFloats.isnormalized`.

The major engineering conclusion is now different from the earlier x4 status:
**the x4 arithmetic kernel is no longer the primary blocker for native linear
algebra.** The current five-pass QW route is both canonical under the audited
corpus and materially faster than upstream mul+add. The next bottleneck is
matrix-kernel data movement, explicit SIMD tiling, packing and threading.

## arXiv:2607.11391 x4 audit correction

The former project diagnosis conflated two different FastTwoSum contracts.

The old diagnostic checked the stronger textbook sufficient condition
`|a| >= |b|` and described violations at the first QW FastTwoSum as an invalid
paper assumption. FPANVerifier and the paper instead use an exponent-order
condition. Reproduction showed:

- former two-pass QW: FastTwoSum gates remain valid under the paper contract,
  but full output non-overlap is not guaranteed;
- current five-pass QW: adds three additional normalization cascades and restores
  a canonical four-limb output in the tested ordinary, wide-exponent and deep
  cancellation corpora;
- 5,000 destructive-cancellation cases: former two-pass had 1,054
  non-normalized outputs; five-pass had 0;
- current five-pass output was bitwise equal to the conservative project x4
  result on the audited corpus;
- all ten concrete five-pass FastTwoSum exponent-contract checks had zero
  violations in the corpus.

The full source-mirrored current-QW FPAN model is expensive enough to exceed a
1,500-second hosted-runner budget when run serially. That is recorded as a proof
budget limit, not a mathematical refutation. The former two-pass expected
refutation pattern was reproduced successfully.

## x4 performance decision

The previous conservative generic `renormalize` finalizer is no longer the
recommended x4 hot path.

After adding the current five-pass QW finalizer and a width-specific x4 public
dispatch, production `fma_fast(Float64x4)` reaches the same performance class as
the standalone paper reproduction.

Representative Ice Lake hosted snapshot, 20,000 scalar FMAs:

| Route | Time |
|---|---:|
| former two-pass QW | ~0.149 ms |
| current five-pass QW | ~0.176 ms |
| production project x4 | ~0.174 ms |
| upstream MultiFloats `x*y+c` | ~0.362 ms |

Thus the current production candidate is roughly **2x faster** than the
separately-rounded upstream scalar composition on that runner while retaining
canonical output in the audited cancellation corpus.

The exact factor varies by hosted CPU, but both Ice Lake and Zen 3 measurements
show the same qualitative result after specialization: current QW is a useful
scalar and SIMD kernel.

## M7 dense linear-algebra baseline

`MultiFloatArithmetic.MFLinearAlgebra` currently provides:

- `mfdot`
- `axpy!`
- `gemv` / `gemv!`
- `gemm` / `gemm!`
- `syrk` / `syrk!`
- `trsv!` / `trsm!`
- `potrf!`

Supported hot types are Float32/Float64 MultiFloat N=2:4.

### Kernel contract

GEMM uses column-major loop order `j -> k -> i`:

```julia
for j in 1:n
    for k in 1:K
        b = B[k,j]
        @simd for i in 1:m
            C[i,j] = fma_fast(A[i,k], b, C[i,j])
        end
    end
end
```

The `i` loop may vectorize because different output rows are independent. One
output entry's reduction over `k` remains ordered. `mfdot` is deliberately not
marked SIMD because its accumulator has a loop-carried dependency.

TRSM and Cholesky use FMA residual updates plus the reproducible
`MultiFloats.div_r` / `sqrt_r` operations.

### Correctness status

Permanent matrix-level validation uses independent 1024-bit BigFloat references.
Linux Julia 1.10/current and current macOS have passed:

- DOT and AXPY;
- GEMV and GEMM, including nontrivial alpha/beta scaling;
- SYRK symmetry/reference accuracy;
- TRSV/TRSM solution accuracy;
- Cholesky reconstruction;
- normalization checks for all generated MultiFloat outputs.

### First performance baseline

Ice Lake, Julia 1.10.11, MultiFloats 3.2.6; comparator is Julia/LinearAlgebra's
generic MultiFloat path.

Float64x2 GEMM:

| n | direct-FMA | generic | speedup |
|---:|---:|---:|---:|
| 16 | 0.007 ms | 0.029 ms | 4.41x |
| 32 | 0.040 ms | 0.214 ms | 5.31x |
| 48 | 0.127 ms | 0.727 ms | 5.71x |

Float64x4 GEMM:

| n | direct-FMA | generic | speedup |
|---:|---:|---:|---:|
| 16 | 0.038 ms | 0.185 ms | 4.92x |
| 32 | 0.274 ms | 1.498 ms | 5.47x |
| 48 | 0.918 ms | 5.017 ms | 5.47x |

GEMV speedups in the first snapshot were about 1.2x for x2 and 2.0x for x4.
The first DOT timing was measured before removal of an unsafe reduction `@simd`
annotation and is therefore not used as an accepted performance result until the
rerun completes.

These results establish a useful baseline before explicit packing or
`MultiFloatVec` microkernels.

## Current comparison with released MultiFloats v3.2.6

The project pins MultiFloats v3.2.6.

Upstream add/multiply remain useful primitives for x2-x4. This project should
not duplicate them without a concrete benefit. The project's current advantages
are instead:

- a direct FMA operation with audited representation behavior;
- faster FMA-based small/medium dense matrix kernels;
- exact higher-width rejection oracles and correctness baselines;
- reproducible proof/corpus infrastructure.

The linear-algebra benchmark shows that using the direct FMA at the matrix level
can produce much larger end-to-end gains than judging scalar add/mul operations
in isolation.

## Higher-width x5-x8 status

x5-x8 remain a separate correctness track.

Current safe multiplication builds an intentionally over-complete expansion, so
x8 is far slower than BigFloat512. This should not be micro-optimized into the
linear-algebra backend.

The higher-width path remains:

1. keep exact/safe implementations as rejection oracles;
2. derive smaller diagonal/truncated fixed-cost multiplication networks;
3. build matching direct-FMA networks without rounded mul-then-add composition;
4. target `MultiFloatVec` and matrix workloads from the start;
5. promote x5/x6 survivors first;
6. lift to x8 only if it becomes competitive with BigFloat512 in intended
   workloads.

## Current decisions

| Component | Decision |
|---|---|
| x2 direct FMA | keep; hot scalar/SIMD primitive |
| x3 direct FMA | keep; proof-driven repaired network |
| x4 former two-pass QW | retain only as audit/performance history |
| x4 current five-pass QW | preferred hot FMA candidate |
| generic x4 renormalize finalizer | no longer preferred hot path |
| x5-x8 `reference_*` | authoritative Experimental rejection oracle |
| x5-x8 safe add/mul/FMA | correctness/reference baselines |
| current x8 safe mul/FMA | never use as current M7 hot path |
| FMA-native DOT/GEMV/GEMM/SYRK | accepted first M7 baseline |
| TRSV/TRSM/POTRF | accepted correctness-first structured baseline |
| explicit packed SIMD GEMM | next priority |

## Immediate next performance milestone

Do not return to generic x4 renormalization tuning.

Next work should optimize the matrix backend while preserving the now-frozen
arithmetic order:

1. explicit `MultiFloatVec` W=4/W=8 row microkernels;
2. B-panel packing and benchmark-derived MR/NR/KC/MC/NC choices;
3. blocked SYRK/TRSM/POTRF;
4. threading over independent output panels, never silently across one ordered
   accumulator reduction;
5. transpose routes and structured/KKT kernels;
6. downstream SDPX.jl A/B on correctness, iterations, wall time, RSS, route and
   certificate residuals.
