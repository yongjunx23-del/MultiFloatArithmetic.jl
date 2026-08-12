module MFLinearAlgebra

using LinearAlgebra
using MultiFloats
using ..MultiFloatArithmetic: fma_fast

import MultiFloats: MultiFloat, div_r, sqrt_r

export axpy!, gemm, gemm!, gemv, gemv!, mfdot, potrf!, syrk, syrk!, trsm!, trsv!

const _SupportedBase = Union{Float32,Float64}

@inline function _check_width(::Type{MultiFloat{T,N}}) where {T,N}
    T <: _SupportedBase || throw(ArgumentError("MFLinearAlgebra supports Float32/Float64 limbs"))
    N in (2, 3, 4) || throw(ArgumentError(
        "MFLinearAlgebra currently supports MultiFloat widths N = 2, 3, or 4"))
    return nothing
end

@inline _coerce(::Type{M}, x::M) where {M} = x
@inline _coerce(::Type{M}, x) where {M} = M(x)

@inline function _scale(x::M, a::M) where {M<:MultiFloat}
    iszero(a) && return zero(M)
    isone(a) && return x
    return fma_fast(a, x, zero(M))
end

@inline function _scale_or_zero(x::M, a::M) where {M<:MultiFloat}
    iszero(a) ? zero(M) : _scale(x, a)
end

"""
    mfdot(x, y)

Compute a MultiFloat dot product with a chained direct-FMA accumulator.

For x2/x3/x4 the accumulator returned by every `fma_fast` step is canonical,
so the next multiply-add starts from a valid normalized expansion. This is the
primitive used by the higher-level matrix kernels in this module.
"""
function mfdot(
    x::AbstractVector{MultiFloat{T,N}},
    y::AbstractVector{MultiFloat{T,N}},
) where {T,N}
    M = MultiFloat{T,N}
    _check_width(M)
    length(x) == length(y) || throw(DimensionMismatch(
        "dot product lengths differ: $(length(x)) and $(length(y))"))

    s = zero(M)
    @inbounds @simd for idx in eachindex(x, y)
        s = fma_fast(x[idx], y[idx], s)
    end
    return s
end

"""
    axpy!(α, x, y)

In-place `y := α*x + y` using one direct MultiFloat FMA per element.
"""
function axpy!(
    α,
    x::AbstractVector{MultiFloat{T,N}},
    y::AbstractVector{MultiFloat{T,N}},
) where {T,N}
    M = MultiFloat{T,N}
    _check_width(M)
    length(x) == length(y) || throw(DimensionMismatch(
        "AXPY lengths differ: $(length(x)) and $(length(y))"))
    a = _coerce(M, α)

    @inbounds @simd for idx in eachindex(x, y)
        y[idx] = fma_fast(a, x[idx], y[idx])
    end
    return y
end

"""
    gemv!(y, A, x; α=1, β=0)

Compute `y := α*A*x + β*y` for dense column-major MultiFloat arrays.

The hot path streams columns of `A`: `k -> i`. With `α == 1`, every matrix
entry contributes through exactly one `fma_fast`, avoiding a separate multiply
and add while preserving a canonical accumulator after each update.
"""
function gemv!(
    y::AbstractVector{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    x::AbstractVector{MultiFloat{T,N}};
    α=one(MultiFloat{T,N}),
    β=zero(MultiFloat{T,N}),
) where {T,N}
    M = MultiFloat{T,N}
    _check_width(M)
    m, k = size(A)
    length(x) == k || throw(DimensionMismatch(
        "A has $k columns but x has length $(length(x))"))
    length(y) == m || throw(DimensionMismatch(
        "A has $m rows but y has length $(length(y))"))
    Base.mightalias(y, x) && throw(ArgumentError("gemv! does not permit y/x aliasing"))

    a = _coerce(M, α)
    b = _coerce(M, β)

    if iszero(b)
        fill!(y, zero(M))
    elseif !isone(b)
        @inbounds @simd for i in eachindex(y)
            y[i] = _scale(y[i], b)
        end
    end

    iszero(a) && return y

    if isone(a)
        @inbounds for p in 1:k
            xp = x[p]
            @simd for i in 1:m
                y[i] = fma_fast(A[i,p], xp, y[i])
            end
        end
    else
        @inbounds for p in 1:k
            xp = _scale(x[p], a)
            @simd for i in 1:m
                y[i] = fma_fast(A[i,p], xp, y[i])
            end
        end
    end
    return y
end

function gemv(
    A::AbstractMatrix{MultiFloat{T,N}},
    x::AbstractVector{MultiFloat{T,N}};
    α=one(MultiFloat{T,N}),
) where {T,N}
    y = zeros(MultiFloat{T,N}, size(A,1))
    return gemv!(y, A, x; α=α, β=zero(MultiFloat{T,N}))
end

"""
    gemm!(C, A, B; α=1, β=0)

Compute `C := α*A*B + β*C` using a streaming column-major direct-FMA kernel.

The default `α=1` route uses loop order `j -> k -> i`: `B[k,j]` is reused
while columns of `A` and `C` are traversed contiguously. This layout is chosen
so LLVM can SIMD-vectorize independent row updates while each scalar MultiFloat
accumulator remains a canonical x2/x3/x4 expansion.
"""
function gemm!(
    C::AbstractMatrix{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}};
    α=one(MultiFloat{T,N}),
    β=zero(MultiFloat{T,N}),
) where {T,N}
    M = MultiFloat{T,N}
    _check_width(M)
    m, k = size(A)
    kb, n = size(B)
    k == kb || throw(DimensionMismatch(
        "inner dimensions differ: size(A,2)=$k and size(B,1)=$kb"))
    size(C) == (m,n) || throw(DimensionMismatch(
        "C has size $(size(C)); expected ($m, $n)"))
    (Base.mightalias(C, A) || Base.mightalias(C, B)) && throw(ArgumentError(
        "gemm! does not permit C to alias A or B"))

    a = _coerce(M, α)
    b = _coerce(M, β)

    if iszero(b)
        fill!(C, zero(M))
    elseif !isone(b)
        @inbounds for j in 1:n
            @simd for i in 1:m
                C[i,j] = _scale(C[i,j], b)
            end
        end
    end

    iszero(a) && return C

    if isone(a)
        @inbounds for j in 1:n
            for p in 1:k
                bpj = B[p,j]
                @simd for i in 1:m
                    C[i,j] = fma_fast(A[i,p], bpj, C[i,j])
                end
            end
        end
    else
        @inbounds for j in 1:n
            for p in 1:k
                bpj = _scale(B[p,j], a)
                @simd for i in 1:m
                    C[i,j] = fma_fast(A[i,p], bpj, C[i,j])
                end
            end
        end
    end
    return C
end

function gemm(
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}};
    α=one(MultiFloat{T,N}),
) where {T,N}
    C = zeros(MultiFloat{T,N}, size(A,1), size(B,2))
    return gemm!(C, A, B; α=α, β=zero(MultiFloat{T,N}))
end

"""
    syrk!(C, A; α=1, β=0, uplo=:L, mirror=true)

Compute the real symmetric rank-k update `C := α*A*A' + β*C`.
Only the requested triangle is accumulated; `mirror=true` copies it to the
opposite triangle after the update.
"""
function syrk!(
    C::AbstractMatrix{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}};
    α=one(MultiFloat{T,N}),
    β=zero(MultiFloat{T,N}),
    uplo::Symbol=:L,
    mirror::Bool=true,
) where {T,N}
    M = MultiFloat{T,N}
    _check_width(M)
    m, k = size(A)
    size(C) == (m,m) || throw(DimensionMismatch(
        "SYRK requires C to be $m x $m; got $(size(C))"))
    uplo in (:L, :U) || throw(ArgumentError("uplo must be :L or :U"))
    Base.mightalias(C, A) && throw(ArgumentError("syrk! does not permit C/A aliasing"))

    a = _coerce(M, α)
    b = _coerce(M, β)

    if uplo === :L
        @inbounds for j in 1:m
            for i in j:m
                C[i,j] = _scale_or_zero(C[i,j], b)
            end
        end
    else
        @inbounds for j in 1:m
            for i in 1:j
                C[i,j] = _scale_or_zero(C[i,j], b)
            end
        end
    end

    if !iszero(a)
        if uplo === :L
            @inbounds for p in 1:k
                for j in 1:m
                    ajp = isone(a) ? A[j,p] : _scale(A[j,p], a)
                    @simd for i in j:m
                        C[i,j] = fma_fast(A[i,p], ajp, C[i,j])
                    end
                end
            end
        else
            @inbounds for p in 1:k
                for j in 1:m
                    ajp = isone(a) ? A[j,p] : _scale(A[j,p], a)
                    @simd for i in 1:j
                        C[i,j] = fma_fast(A[i,p], ajp, C[i,j])
                    end
                end
            end
        end
    end

    if mirror
        if uplo === :L
            @inbounds for j in 1:m, i in j+1:m
                C[j,i] = C[i,j]
            end
        else
            @inbounds for j in 1:m, i in 1:j-1
                C[j,i] = C[i,j]
            end
        end
    end
    return C
end

function syrk(
    A::AbstractMatrix{MultiFloat{T,N}};
    α=one(MultiFloat{T,N}),
    uplo::Symbol=:L,
) where {T,N}
    C = zeros(MultiFloat{T,N}, size(A,1), size(A,1))
    return syrk!(C, A; α=α, β=zero(MultiFloat{T,N}), uplo=uplo, mirror=true)
end

"""
    trsv!(A, x; uplo=:L, unitdiag=false)

Solve `A*x = b` in place for a triangular MultiFloat matrix `A`, overwriting
`x` with the solution. Division uses `MultiFloats.div_r` for reproducibility.
"""
function trsv!(
    A::AbstractMatrix{MultiFloat{T,N}},
    x::AbstractVector{MultiFloat{T,N}};
    uplo::Symbol=:L,
    unitdiag::Bool=false,
) where {T,N}
    M = MultiFloat{T,N}
    _check_width(M)
    n, n2 = size(A)
    n == n2 || throw(DimensionMismatch("triangular solve requires a square A"))
    length(x) == n || throw(DimensionMismatch(
        "A is $n x $n but x has length $(length(x))"))
    uplo in (:L, :U) || throw(ArgumentError("uplo must be :L or :U"))

    if uplo === :L
        @inbounds for i in 1:n
            s = x[i]
            for p in 1:i-1
                s = fma_fast(-A[i,p], x[p], s)
            end
            x[i] = unitdiag ? s : div_r(s, A[i,i])
        end
    else
        @inbounds for i in n:-1:1
            s = x[i]
            for p in i+1:n
                s = fma_fast(-A[i,p], x[p], s)
            end
            x[i] = unitdiag ? s : div_r(s, A[i,i])
        end
    end
    return x
end

"""
    trsm!(A, B; uplo=:L, unitdiag=false)

Solve `A*X = B` in place, overwriting `B` with `X`.
"""
function trsm!(
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}};
    uplo::Symbol=:L,
    unitdiag::Bool=false,
) where {T,N}
    M = MultiFloat{T,N}
    _check_width(M)
    n, n2 = size(A)
    n == n2 || throw(DimensionMismatch("triangular solve requires a square A"))
    size(B,1) == n || throw(DimensionMismatch(
        "A is $n x $n but B has $(size(B,1)) rows"))
    uplo in (:L, :U) || throw(ArgumentError("uplo must be :L or :U"))
    Base.mightalias(A, B) && throw(ArgumentError("trsm! does not permit A/B aliasing"))

    nrhs = size(B,2)
    if uplo === :L
        @inbounds for j in 1:nrhs
            for i in 1:n
                s = B[i,j]
                for p in 1:i-1
                    s = fma_fast(-A[i,p], B[p,j], s)
                end
                B[i,j] = unitdiag ? s : div_r(s, A[i,i])
            end
        end
    else
        @inbounds for j in 1:nrhs
            for i in n:-1:1
                s = B[i,j]
                for p in i+1:n
                    s = fma_fast(-A[i,p], B[p,j], s)
                end
                B[i,j] = unitdiag ? s : div_r(s, A[i,i])
            end
        end
    end
    return B
end

"""
    potrf!(A; uplo=:L, clean=true)

In-place Cholesky factorization of a real symmetric positive-definite
MultiFloat matrix. The default computes a lower-triangular `L` such that
`A_original ≈ L*L'`. Square root and division use the reproducible MultiFloats
variants `sqrt_r` and `div_r`.

The upper-triangular route factors by symmetry and stores `U` with
`A_original ≈ U'*U`.
"""
function potrf!(
    A::AbstractMatrix{MultiFloat{T,N}};
    uplo::Symbol=:L,
    clean::Bool=true,
) where {T,N}
    M = MultiFloat{T,N}
    _check_width(M)
    n, n2 = size(A)
    n == n2 || throw(DimensionMismatch("Cholesky factorization requires a square matrix"))
    uplo in (:L, :U) || throw(ArgumentError("uplo must be :L or :U"))

    if uplo === :L
        @inbounds for j in 1:n
            d = A[j,j]
            for p in 1:j-1
                d = fma_fast(-A[j,p], A[j,p], d)
            end
            A[j,j] = sqrt_r(d)

            for i in j+1:n
                s = A[i,j]
                for p in 1:j-1
                    s = fma_fast(-A[i,p], A[j,p], s)
                end
                A[i,j] = div_r(s, A[j,j])
            end
        end
        if clean
            @inbounds for j in 1:n, i in 1:j-1
                A[i,j] = zero(M)
            end
        end
    else
        @inbounds for j in 1:n
            d = A[j,j]
            for p in 1:j-1
                d = fma_fast(-A[p,j], A[p,j], d)
            end
            A[j,j] = sqrt_r(d)

            for k in j+1:n
                s = A[j,k]
                for p in 1:j-1
                    s = fma_fast(-A[p,j], A[p,k], s)
                end
                A[j,k] = div_r(s, A[j,j])
            end
        end
        if clean
            @inbounds for j in 1:n, i in j+1:n
                A[i,j] = zero(M)
            end
        end
    end
    return A
end

end # module MFLinearAlgebra
