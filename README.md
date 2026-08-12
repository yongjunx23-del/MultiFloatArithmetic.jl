# MultiFloatArithmetic.jl

[![CI](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/yongjunx23-del/MultiFloatArithmetic.jl/actions/workflows/ci.yml)

Verification-oriented arithmetic and linear-algebra research kernels for
[`MultiFloats.jl`](https://github.com/dzhang314/MultiFloats.jl).

The project develops fixed-length MultiFloat arithmetic together with native
DOT/GEMV/GEMM/SYRK/TRSM/factorization kernels for downstream high-precision
solvers such as SDPX.jl. Correctness and representation contracts come before
instruction-count wins.

## Current stage

The project now has two distinct tiers.

### Hot x2-x4 arithmetic and linear algebra

- x2/x3 use fixed-cost direct-FMA networks with structural verification;
- x4 uses the audited five-pass QW normalization from the current
  arXiv:2607.11391 algorithm rather than the former generic `renormalize`
  fallback;
- the x4 public dispatch is width-specialized so Julia emits the same class of
  code as the standalone QW reproduction;
- `MFLinearAlgebra` provides FMA-native `mfdot`, AXPY, GEMV, GEMM, SYRK,
  TRSV/TRSM and Cholesky for Float32/Float64 N=2:4;
- chained linear-algebra accumulators are checked against 1024-bit BigFloat and
  must remain `MultiFloats.isnormalized`.

On the first Ice Lake linear-algebra snapshot, the simple column-major FMA GEMM
was about **5.3-5.7x faster** than Julia/LinearAlgebra's generic MultiFloat path
for Float64x2 and about **4.9-5.5x faster** for Float64x4 at n=16:48, before any
panel packing or hand-written `MultiFloatVec` microkernel.

### x5-x8 correctness baselines

- exact-rational Experimental add/sub/mul/FMA oracles cover x5-x8;
- `add_safe` and `mul_safe` provide no-FastTwoSum correctness families;
- `fma_safe` forms one exact direct `x*y+c` expansion before N-limb truncation;
- the higher-limb implementations are correctness references, not hot kernels.

In particular, the intentionally over-complete x8 multiplication is much slower
than BigFloat512 and must be replaced by a fixed-cost network before x5-x8 enter
the linear-algebra backend.

See [STATUS.md](STATUS.md), [docs/LINEAR_ALGEBRA.md](docs/LINEAR_ALGEBRA.md),
[docs/NUMERICAL_CONTRACT.md](docs/NUMERICAL_CONTRACT.md),
[docs/REFERENCE_X5_X8.md](docs/REFERENCE_X5_X8.md),
[docs/ADD_SAFE_X5_X8.md](docs/ADD_SAFE_X5_X8.md),
[docs/MUL_SAFE_X5_X8.md](docs/MUL_SAFE_X5_X8.md),
[docs/FMA_SAFE_X5_X8.md](docs/FMA_SAFE_X5_X8.md), and
[docs/BIGFLOAT_MULTIFLOAT_COMPARISON.md](docs/BIGFLOAT_MULTIFLOAT_COMPARISON.md).

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

## MultiFloat-native linear algebra

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

The first dense backend supports Float32/Float64 x2-x4. GEMV/GEMM/SYRK use
ordered chained FMA accumulators, while TRSM/Cholesky additionally use
`MultiFloats.div_r` and `sqrt_r` for reproducible division and square root.

No methods are added to external `LinearAlgebra.dot`/`mul!`; the namespaced API
avoids type piracy. See [docs/LINEAR_ALGEBRA.md](docs/LINEAR_ALGEBRA.md) for the
loop contract, validation, benchmarks and next microkernel layer.

## Higher-limb exact oracle

```julia
const E = MultiFloatArithmetic.Experimental
T8 = MultiFloat{Float64,8}
x8 = T8(BigFloat("1.2345678901234567890123456789"))
y8 = T8(BigFloat("0.9876543210987654321098765432"))
c8 = T8(BigFloat("-0.1250000000000000000000000001"))

oracle_sum = E.reference_add(x8, y8)
oracle_product = E.reference_mul(x8, y8)
oracle_fma = E.reference_fma(x8, y8, c8)
```

`reference_add/sub/mul/fma` support Float32/Float64 N=5:8 and deliberately use
exact `Rational{BigInt}` arithmetic before one final pack.

## Safe higher-limb add/mul

```julia
safe_sum = E.add_safe(x8, y8)
safe_product = E.mul_safe(x8, y8)
```

Both families preserve and normalize over-complete exact expansions before
N-limb truncation and introduce no FastTwoSum assumptions. They are correctness
baselines, not hot APIs.

## Direct Float64x5-x8 FMA baseline

```julia
safe_fma = E.fma_safe(x8, y8, c8)
# Width-specific wrappers E.fma5_safe ... E.fma8_safe are also available.
```

For width N, `fma_safe` keeps all `2N^2` exact TwoProd product/residual components
and all N limbs of `c`, canonicalizes the resulting `2N^2+N` component direct
`x*y+c` expansion, fully renormalizes it, and only then truncates to N limbs.
It does not call the rounded N-limb `mul_safe` result and does not use FastTwoSum.

The permanent corpus requires exact full-expansion head+tail accounting,
normalized output, bitwise x/y symmetry, and bitwise equality with
`reference_fma`, including deep destructive cancellation.

First maximum observed `|err|/(u^N(|xy|+|c|))` values were approximately:

- x5: 0.04928
- x6: 0.02276
- x7: 0.00452
- x8: 0.000826

The rounded composition `add_safe(mul_safe(x,y),c)` disagreed with the direct
exact-FMA oracle in every destructive-cancellation case tested at x5-x8, while
the direct path had zero oracle mismatches. The current high-limb direct FMA is
therefore a correctness baseline that still needs a smaller fixed-cost network.

## BigFloat versus MultiFloat comparison

A permanent benchmark compares:

- BigFloat128 vs Float64x2 (~106 nominal significand bits),
- BigFloat256 vs Float64x4 (~212 bits),
- BigFloat512 vs Float64x8 (~424 bits).

The BigFloat settings have about 20.8% more nominal precision, so speed and
precision are reported together. Original upstream x2/x4 arithmetic has large
throughput advantages over BigFloat128/256, whereas the current correctness-only
x8 multiplication is far slower than BigFloat512.

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
include("benchmark/paper2607_v4_audit.jl")
include("benchmark/linear_algebra.jl")
include("benchmark/add5_safe.jl")
include("benchmark/add6_add8_safe.jl")
include("benchmark/mul5_safe.jl")
include("benchmark/mul6_mul8_safe.jl")
include("benchmark/fma5_safe.jl")
include("benchmark/fma6_fma8_safe.jl")
include("benchmark/bigfloat_multifloat_comparison.jl")
include("benchmark/division_digits.jl")
include("benchmark/codegen.jl")
```

## Roadmap

1. Add explicit `MultiFloatVec` GEMM microkernels (W=4/8 first), B-panel packing,
   and benchmark-derived MC/KC/NC blocking while preserving one accumulator's
   ordered reduction.
2. Build blocked SYRK/TRSM/POTRF from the same kernels and add threaded outer
   panel scheduling.
3. Add transpose routes, LDLT/pivoting, solve/refine primitives and KKT-oriented
   structured kernels for optimization workloads.
4. Derive smaller formally justified x5-x8 arithmetic networks; only competitive
   survivors enter the linear-algebra layer.
5. A/B the backend inside SDPX.jl on correctness, iterations, time, RSS and
   certification residuals.

Research software; passing CI is reproducible evidence, not a production/formal
guarantee.

MIT.
