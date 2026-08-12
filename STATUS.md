# Project status

Status date: **2026-08-12**  
Current milestone: **M7.2 — explicit-SIMD GEMM and blocked Float64x4 Cholesky**

## Executive assessment

The project is now a usable fixed-width MultiFloat arithmetic plus dense
linear-algebra research layer rather than an arithmetic-only prototype.

Accepted direction:

- x2/x3 fixed-cost direct-FMA networks with pinned structural verification;
- x4 current five-pass QW direct FMA, reproduced from arXiv:2607.11391 and
  width-specialized for Julia code generation;
- x5-x8 exact/safe arithmetic remains Experimental and is not used in M7 hot
  paths;
- `MFLinearAlgebra` provides FMA-native DOT/AXPY/GEMV/GEMM/SYRK/TRSV/TRSM and
  Cholesky for Float32/Float64 x2-x4;
- dense Float64x4 GEMM now has an explicit `MultiFloatVec` MR=8, NR=2 route;
- dense lower Float64x4 Cholesky now has a bitwise-preserving Vec8 blocked route
  for `n >= 32` with block size 8;
- matrix outputs remain checked against high-precision references and must remain
  `MultiFloats.isnormalized`.

The arithmetic layer is no longer the main x4 performance blocker. Current M7
work is primarily about data reuse, vector width, structured blocking and, later,
threading while preserving each ordered chained-FMA reduction.

## arXiv:2607.11391 x4 audit

The former project diagnosis conflated two different FastTwoSum contracts. The
old diagnostic used the stronger textbook sufficient condition `|a| >= |b|`;
the paper and FPANVerifier use an exponent-order condition.

Reproduction established:

- former 146-flop/two-pass QW: FastTwoSum gates remain valid under the paper
  contract, but full output non-overlap is not guaranteed;
- current 176-flop/five-pass QW: all ten concrete FastTwoSum exponent-contract
  checks had zero violations in the audit corpus;
- 5,000 destructive-cancellation cases: the two-pass result was non-normalized
  in 1,054 cases; the five-pass result had zero normalization failures;
- current five-pass output was bitwise equal to the conservative project x4
  result on ordinary, wide-exponent and cancellation audit corpora.

The full current-QW universal FPAN proof is expensive enough to exceed the
hosted-runner budget when run serially. That remains a proof-budget limitation,
not a mathematical refutation.

## x4 arithmetic performance

The five-pass QW route is the production x4 FMA. A representative Zen 3 run,
20,000 scalar operations:

| Route | Time |
|---|---:|
| former two-pass QW | 0.165 ms |
| current five-pass QW / project `fma_fast` | 0.218-0.219 ms |
| upstream MultiFloats `x*y+c` | 0.373-0.374 ms |

Thus the current canonical direct FMA is about **1.7x faster** than separately
rounded upstream `x*y+c` on that runner. Vec2/Vec4/Vec8 direct FMA remained about
1.9x/1.9x/1.5x faster than upstream in the same audit run.

## M7 dense linear algebra

`MultiFloatArithmetic.MFLinearAlgebra` provides:

- `mfdot`
- `axpy!`
- `gemv` / `gemv!`
- `gemm` / `gemm!`
- `syrk` / `syrk!`
- `trsv!` / `trsm!`
- `potrf!`

Supported hot types are Float32/Float64 MultiFloat N=2:4.

### Numerical scheduling contract

A single output accumulator retains its chronological reduction order. SIMD is
used only across independent rows/outputs. The project does not silently reassociate
one dot-product/FMA chain to obtain speed.

`mfdot` is therefore intentionally ordered. GEMM and structured trailing updates
vectorize independent rows. TRSM/POTRF continue to use reproducible `div_r` and
`sqrt_r` at their scalar dependency points.

### Matrix correctness status

Permanent matrix tests use independent 1024-bit BigFloat references and run on
Linux Julia 1.10/current and current macOS. They cover:

- DOT/AXPY;
- GEMV/GEMM, including nontrivial alpha/beta routes;
- SYRK symmetry/reference accuracy;
- TRSV/TRSM solution accuracy;
- Cholesky reconstruction;
- normalization of all generated MultiFloat outputs.

Optimized dense routes carry an additional stronger regression gate: where the
optimization only changes scheduling/vector grouping, the optimized output must
be **bitwise `===`** to its accepted scalar/streaming baseline.

## M7.1 explicit-SIMD GEMM

The first GEMM baseline used `j -> k -> i` streaming and compiler `@simd`. A
permanent A/B then compared explicit `MultiFloatVec` widths and two-column reuse.

Decision:

- explicit Vec4: rejected;
- explicit x2 route: not promoted because compiler streaming was already very
  strong;
- Float64x4 MR=8, NR=2: promoted for dense `Matrix`, unit alpha, `m >= 8`,
  `n >= 2`, `k >= 4`;
- all other shapes/types/scaling routes retain the streaming fallback.

Zen 3 public Float64x4 GEMM after promotion:

| n | public | streaming | public speedup | generic `mul!` | public/generic |
|---:|---:|---:|---:|---:|---:|
| 16 | 0.038 ms | 0.045 ms | 1.19x | 0.290 ms | 7.66x |
| 32 | 0.277 ms | 0.347 ms | 1.25x | 2.351 ms | 8.48x |
| 48 | 0.923 ms | 1.188 ms | 1.29x | 7.890 ms | 8.55x |
| 64 | 2.205 ms | 2.796 ms | 1.27x | 18.653 ms | 8.46x |

The explicit route is elementwise bitwise equal to the ordered streaming kernel.

On the same run the existing x2 streaming GEMM was roughly 7.3-9.3x faster than
generic `mul!`; explicit x2 vector packing therefore remains unnecessary.

## M7.2 blocked Float64x4 Cholesky

The lower Cholesky kernel was reorganized into diagonal panels and trailing Schur
updates. The optimization has a deliberately stronger-than-usual numerical
contract:

- every destination element sees the same `p = 1,2,...` FMA chronology;
- FMA operand roles are preserved exactly;
- Vec8 only groups independent rows;
- `sqrt_r` and `div_r` operations are unchanged;
- the full output matrix must be bitwise equal to the unblocked baseline.

Block sizes 8/16/24/32 were tested. BS=8 was consistently best and is promoted
only for `Matrix{Float64x4}`, lower factorization, `n >= 32`; small matrices,
upper factorization and all other types use the unblocked fallback.

Promotion-run Zen 3 results:

| n | public BS=8 | unblocked | speedup |
|---:|---:|---:|---:|
| 32 | 0.212 ms | 0.365 ms | **1.72x** |
| 64 | 1.058 ms | 2.468 ms | **2.33x** |
| 96 | 2.850 ms | 7.801 ms | **2.74x** |

Permanent n=17/33 tests, a benchmark n=73 remainder case and all measured block
sizes were whole-matrix bitwise equal to the unblocked result. Linux 1.10/current
and macOS correctness jobs are green with the public dispatch enabled.

## Current comparison with MultiFloats v3.2.6

The project pins MultiFloats v3.2.6. Upstream add/multiply remain useful x2-x4
primitives and are not duplicated without evidence of benefit. This project now
adds value primarily through:

- audited direct FMA;
- explicit fixed-width SIMD scheduling;
- substantially faster dense GEMM and structured Cholesky paths;
- exact higher-width rejection oracles/correctness baselines;
- reproducible proof, correctness and benchmark infrastructure.

The matrix-level results show why scalar add/mul microbenchmarks alone are not a
sufficient measure of the fixed-width MultiFloat design.

## Higher-width x5-x8 status

x5-x8 remain a separate correctness track. Current safe multiplication builds an
intentionally over-complete expansion, and x8 is far slower than BigFloat512.
It must not be inserted into the M7 backend as-is.

The higher-width path remains:

1. keep exact/safe implementations as rejection oracles;
2. derive smaller diagonal/truncated multiplication networks;
3. build matching direct-FMA networks without rounded mul-then-add composition;
4. target `MultiFloatVec` and matrix workloads from the start;
5. promote x5/x6 survivors first;
6. lift to x8 only if intended workloads become competitive with BigFloat512.

## Current decisions

| Component | Decision |
|---|---|
| x2 direct FMA | keep; hot scalar/SIMD primitive |
| x3 direct FMA | keep; proof-driven repaired network |
| x4 former two-pass QW | audit/performance history only |
| x4 current five-pass QW | accepted hot direct FMA |
| x5-x8 `reference_*` | authoritative Experimental rejection oracle |
| x5-x8 safe add/mul/FMA | correctness/reference baselines |
| current x8 safe mul/FMA | never use as M7 hot path |
| DOT/AXPY/GEMV/SYRK | accepted M7 baseline |
| Float64x4 dense GEMM MR8xNR2 | accepted M7.1 fast path |
| TRSV/TRSM | accepted correctness-first structured baseline |
| Float64x4 lower POTRF BS=8, n>=32 | accepted M7.2 fast path |
| upper/small/other POTRF | unblocked fallback |

## Immediate next performance milestone

M7.3 should optimize structured kernels that can reuse the same proven scheduling
principle:

1. explicit Vec8 lower/upper SYRK A/B with bitwise-equality gates;
2. TRSM optimization across independent RHS columns only after verifying any
   vector-division route is lane-wise equivalent to scalar `div_r`;
3. benchmark whether B/panel packing adds value beyond the current MR8xNR2 GEMM;
4. thread only independent output panels; never parallelize one ordered reduction
   without an explicit new numerical contract;
5. add transpose/structured/KKT kernels needed by downstream SDPX.jl;
6. run downstream solver A/B on correctness, iterations, wall time, RSS, route
   and certificate residuals.
