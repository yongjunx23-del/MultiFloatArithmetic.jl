# FMA-native MultiFloat linear algebra

`MultiFloatArithmetic.MFLinearAlgebra` is the first dense CPU linear-algebra
layer built directly on this project's verified/validated fixed-width FMA
networks.

The design goal is not to wrap generic Julia matrix multiplication. The goal is
to make the MultiFloat arithmetic contract visible to the linear-algebra kernel:
for x2/x3/x4, every chained accumulator is updated with `fma_fast`, and every
update returns a canonical expansion suitable as the next FMA input.

## Supported arithmetic

Current hot linear algebra supports:

- `MultiFloat{Float32,N}` and `MultiFloat{Float64,N}`;
- `N = 2, 3, 4`;
- real dense arrays;
- CPU execution.

x5-x8 are intentionally not routed into these kernels yet. Their current
addition/multiplication/FMA implementations are correctness baselines, not
competitive hot arithmetic.

## API

```julia
using MultiFloats
using MultiFloatArithmetic

const LA = MultiFloatArithmetic.MFLinearAlgebra
T = Float64x4

A = rand(T, 64, 64)
B = rand(T, 64, 64)
C = LA.gemm(A, B)

x = rand(T, 64)
y = LA.gemv(A, x)
d = LA.mfdot(x, y)

S = LA.syrk(A)

L = copy(S)
LA.potrf!(L; uplo=:L)

rhs = rand(T, 64, 4)
LA.trsm!(L, rhs; uplo=:L)
```

Available routines:

- `mfdot(x, y)`
- `axpy!(α, x, y)`
- `gemv(A, x)` / `gemv!(y, A, x; α=1, β=0)`
- `gemm(A, B)` / `gemm!(C, A, B; α=1, β=0)`
- `syrk(A)` / `syrk!(C, A; α=1, β=0, uplo=:L, mirror=true)`
- `trsv!(A, x; uplo=:L, unitdiag=false)`
- `trsm!(A, B; uplo=:L, unitdiag=false)`
- `potrf!(A; uplo=:L, clean=true)`

The module deliberately does **not** add methods to `LinearAlgebra.dot` or
`LinearAlgebra.mul!`. Both those functions and the `MultiFloat` type are owned by
other modules, so doing that here would be type piracy.

## GEMM kernel

The default `α = 1` GEMM hot loop is conceptually:

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

This loop order is intentional for Julia's column-major storage:

1. `B[k,j]` is loaded once and reused down a column;
2. `A[:,k]` and `C[:,j]` are traversed contiguously;
3. different `i` updates are independent and may be SIMD-vectorized;
4. each individual `C[i,j]` accumulation remains an ordered chained FMA over
   `k`.

This is particularly important for x4. The current five-pass QW FMA restores a
canonical four-limb expansion after each update, so the next iteration satisfies
the same normalized-input representation contract. The former two-pass QW
network could lose full non-overlap after cancellation and therefore was not a
safe basis for an unrestricted chained GEMM accumulator.

## Reduction ordering

`mfdot` is deliberately **not** marked `@simd` because it has a loop-carried
MultiFloat accumulator. Its reduction order is part of the numerical result.

The row loops in AXPY/GEMV/GEMM/SYRK *are* marked `@simd`: each lane updates a
different output entry and therefore has no cross-iteration arithmetic
dependency.

## α and β semantics

The `α = 1` routes are the primary hot paths.

For non-unit `α`, GEMV/GEMM scale the reused vector/matrix operand once with a
direct MultiFloat FMA against zero and then use the scaled value in the inner
FMA loop. This avoids multiplying every output entry by `α` separately. As with
BLAS implementations generally, the finite-precision operation order is part of
the implementation and is not a correctly rounded evaluation of the symbolic
real expression.

`β` is applied once to the destination before accumulation.

## Triangular solve and Cholesky

TRSV/TRSM use FMA residual updates:

```julia
s = fma_fast(-A[i,k], x[k], s)
```

and use `MultiFloats.div_r` for the diagonal division.

`potrf!` uses the same FMA residual pattern for diagonal and panel updates and
uses the reproducible `MultiFloats.sqrt_r` / `div_r` primitives. This keeps the
non-FMA operations aligned with MultiFloats' cross-platform reproducible
variants.

The current Cholesky implementation is unblocked. It is a correctness-first
factorization baseline that already exposes the right arithmetic primitives for
a future blocked POTRF built from SYRK/GEMM/TRSM.

## Independent correctness validation

`test/linear_algebra.jl` evaluates x2 and x4 against independent 1024-bit
BigFloat references. The permanent checks cover:

- DOT;
- AXPY;
- GEMV, including nontrivial `α/β`;
- GEMM, including nontrivial `α/β`;
- SYRK symmetry and reference accuracy;
- TRSV/TRSM solution accuracy;
- Cholesky reconstruction;
- `MultiFloats.isnormalized` on all generated results.

The tests run under Julia 1.10/current on Linux and current Julia on macOS.

## First hosted performance snapshot

Ice Lake client, Julia 1.10.11, MultiFloats 3.2.6. The baseline comparator is
Julia/LinearAlgebra's generic MultiFloat path.

### Float64x2

| Kernel | Size | direct-FMA | generic | speedup |
|---|---:|---:|---:|---:|
| DOT | 20,000 | 0.100 ms* | 0.205 ms | 2.06x* |
| GEMV | 256x128 | 0.078 ms | 0.093 ms | 1.19x |
| GEMM | 16 | 0.007 ms | 0.029 ms | 4.41x |
| GEMM | 32 | 0.040 ms | 0.214 ms | 5.31x |
| GEMM | 48 | 0.127 ms | 0.727 ms | 5.71x |

### Float64x4

| Kernel | Size | direct-FMA | generic | speedup |
|---|---:|---:|---:|---:|
| DOT | 20,000 | 0.617 ms* | 1.002 ms | 1.62x* |
| GEMV | 256x128 | 0.291 ms | 0.583 ms | 2.00x |
| GEMM | 16 | 0.038 ms | 0.185 ms | 4.92x |
| GEMM | 32 | 0.274 ms | 1.498 ms | 5.47x |
| GEMM | 48 | 0.918 ms | 5.017 ms | 5.47x |

`*` The first DOT numbers were measured before removal of an unsafe `@simd`
annotation on the ordered reduction. They are retained only as historical
context and must be replaced by the next benchmark run; GEMV/GEMM numbers are
unaffected by that correction.

The important baseline result is already clear: even without packing or a
hand-written SIMD microkernel, the chained-FMA GEMM is roughly five times faster
than the generic MultiFloat matrix path at these small/medium sizes.

## Next performance layer

The current kernel is intentionally simple enough to audit. The next layer
should preserve its numerical contract while changing data movement and SIMD
execution:

1. **Explicit `MultiFloatVec` microkernels** for Float64x2/x4, starting with W=4
   and W=8 row tiles.
2. **B-panel packing** so the inner kernel reuses contiguous packed B data and
   reduces indexing/dispatch overhead.
3. **Blocked GEMM** with small MRxNR register tiles and KC/NC/MC panels chosen
   from benchmark data rather than hard-coded BLAS assumptions.
4. **Blocked SYRK and POTRF** built from the same GEMM/SYRK/TRSM kernels.
5. **Threading over independent output-column/panel tiles**, never over one
   accumulator's reduction dimension unless a separately specified reproducible
   reduction tree is introduced.
6. Transposed GEMM/GEMV routes and structured matrix wrappers.
7. LDLT/pivoted factorization and solve/refinement hooks for optimization/KKT
   workloads.
8. x5-x8 integration only after fixed-cost higher-limb arithmetic becomes
   competitive.

This module is therefore the linear-algebra correctness/performance baseline,
not the final BLAS replacement.
