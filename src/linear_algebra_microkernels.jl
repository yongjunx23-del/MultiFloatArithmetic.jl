# Explicit MultiFloatVec microkernels layered on top of the auditable scalar
# MFLinearAlgebra baseline. Keep the fast path deliberately narrow: current
# evidence supports dense Float64x4 matrices on the hosted Zen3 runner, while
# x2 already auto-vectorizes extremely well and AbstractMatrix layouts may not
# provide contiguous columns.

@eval MFLinearAlgebra begin

import MultiFloats: MultiFloatVec

@inline function _load_mfvec(
    ::Type{MultiFloatVec{W,T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    i0::Int,
    j::Int,
) where {W,T,N}
    return MultiFloatVec{W,T,N}(
        ntuple(lane -> A[i0 + lane - 1, j], Val(W)),
    )
end

@inline function _store_mfvec!(
    A::AbstractMatrix{MultiFloat{T,N}},
    i0::Int,
    j::Int,
    v::MultiFloatVec{W,T,N},
) where {W,T,N}
    @inbounds for lane in 1:W
        A[i0 + lane - 1, j] = v[lane]
    end
    return nothing
end

# Source-equivalent fallback to the public baseline. Keeping it private lets the
# benchmark compare the promoted dense path against the former streaming kernel
# without changing arithmetic semantics or relying on `invoke` with keywords.
function _gemm_streaming_kernel!(
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

"""
    _gemm_vec_accumulate!(C, A, B, Val(W))

Single-column explicit-SIMD diagnostic retained for A/B. It keeps W output rows
in one `MultiFloatVec` accumulator across the complete k reduction. Each lane
sees the same `p = 1:k` reduction order as the scalar streaming kernel.
"""
function _gemm_vec_accumulate!(
    C::AbstractMatrix{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}},
    ::Val{W},
) where {T,N,W}
    M = MultiFloat{T,N}
    V = MultiFloatVec{W,T,N}
    _check_width(M)

    m, k = size(A)
    kb, n = size(B)
    k == kb || throw(DimensionMismatch(
        "inner dimensions differ: size(A,2)=$k and size(B,1)=$kb"))
    size(C) == (m,n) || throw(DimensionMismatch(
        "C has size $(size(C)); expected ($m, $n)"))
    W > 0 || throw(ArgumentError("SIMD width must be positive"))

    full_rows = m - rem(m, W)

    @inbounds for j in 1:n
        for i0 in 1:W:full_rows
            cv = _load_mfvec(V, C, i0, j)
            for p in 1:k
                av = _load_mfvec(V, A, i0, p)
                cv = fma_fast(av, V(B[p,j]), cv)
            end
            _store_mfvec!(C, i0, j, cv)
        end

        for i in full_rows + 1:m
            s = C[i,j]
            for p in 1:k
                s = fma_fast(A[i,p], B[p,j], s)
            end
            C[i,j] = s
        end
    end
    return C
end

"""
    _gemm_vec2col_accumulate!(C, A, B, Val(W))

MR=W, NR=2 register-blocked explicit-SIMD kernel. For every pair of output
columns, the W-row vector from A is loaded once per k step and reused to update
two independent MultiFloatVec accumulators. The arithmetic reduction order for
each output element remains exactly `p = 1:k`, so the result is bitwise
identical to the scalar streaming baseline.

Odd output columns and row remainders use the same reduction order with small
single-column/scalar tails; no padding or changed summation tree is introduced.
"""
function _gemm_vec2col_accumulate!(
    C::AbstractMatrix{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}},
    ::Val{W},
) where {T,N,W}
    M = MultiFloat{T,N}
    V = MultiFloatVec{W,T,N}
    _check_width(M)

    m, k = size(A)
    kb, n = size(B)
    k == kb || throw(DimensionMismatch(
        "inner dimensions differ: size(A,2)=$k and size(B,1)=$kb"))
    size(C) == (m,n) || throw(DimensionMismatch(
        "C has size $(size(C)); expected ($m, $n)"))
    W > 0 || throw(ArgumentError("SIMD width must be positive"))

    full_rows = m - rem(m, W)
    paired_cols = n - rem(n, 2)

    @inbounds for j in 1:2:paired_cols
        for i0 in 1:W:full_rows
            c0 = _load_mfvec(V, C, i0, j)
            c1 = _load_mfvec(V, C, i0, j + 1)
            for p in 1:k
                av = _load_mfvec(V, A, i0, p)
                c0 = fma_fast(av, V(B[p,j]), c0)
                c1 = fma_fast(av, V(B[p,j + 1]), c1)
            end
            _store_mfvec!(C, i0, j, c0)
            _store_mfvec!(C, i0, j + 1, c1)
        end

        for i in full_rows + 1:m
            s0 = C[i,j]
            s1 = C[i,j + 1]
            for p in 1:k
                aip = A[i,p]
                s0 = fma_fast(aip, B[p,j], s0)
                s1 = fma_fast(aip, B[p,j + 1], s1)
            end
            C[i,j] = s0
            C[i,j + 1] = s1
        end
    end

    if isodd(n)
        j = n
        @inbounds for i0 in 1:W:full_rows
            cv = _load_mfvec(V, C, i0, j)
            for p in 1:k
                av = _load_mfvec(V, A, i0, p)
                cv = fma_fast(av, V(B[p,j]), cv)
            end
            _store_mfvec!(C, i0, j, cv)
        end
        @inbounds for i in full_rows + 1:m
            s = C[i,j]
            for p in 1:k
                s = fma_fast(A[i,p], B[p,j], s)
            end
            C[i,j] = s
        end
    end

    return C
end

@inline function _prepare_gemm_vec!(
    C::AbstractMatrix{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}},
    α,
    β,
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
    isone(a) || return (false, a, b)

    if iszero(b)
        fill!(C, zero(M))
    elseif !isone(b)
        @inbounds for j in 1:n
            @simd for i in 1:m
                C[i,j] = _scale(C[i,j], b)
            end
        end
    end
    return (true, a, b)
end

function _gemm_vec!(
    C::AbstractMatrix{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}},
    width::Val{W};
    α=one(MultiFloat{T,N}),
    β=zero(MultiFloat{T,N}),
) where {T,N,W}
    prepared, a, b = _prepare_gemm_vec!(C, A, B, α, β)
    prepared || return _gemm_streaming_kernel!(C, A, B; α=a, β=b)
    return _gemm_vec_accumulate!(C, A, B, width)
end

function _gemm_vec2col!(
    C::AbstractMatrix{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}},
    width::Val{W};
    α=one(MultiFloat{T,N}),
    β=zero(MultiFloat{T,N}),
) where {T,N,W}
    prepared, a, b = _prepare_gemm_vec!(C, A, B, α, β)
    prepared || return _gemm_streaming_kernel!(C, A, B; α=a, β=b)
    return _gemm_vec2col_accumulate!(C, A, B, width)
end

# Promoted dense Float64x4 route. Two independent hosted Zen3 runs measured
# ~15-17% speedup over the already-fast streaming kernel for n=16:64, while
# preserving bitwise equality. Keep the gate conservative so tiny matrices and
# non-unit α stay on the simpler baseline.
function gemm!(
    C::Matrix{MultiFloat{Float64,4}},
    A::Matrix{MultiFloat{Float64,4}},
    B::Matrix{MultiFloat{Float64,4}};
    α=one(MultiFloat{Float64,4}),
    β=zero(MultiFloat{Float64,4}),
)
    M = MultiFloat{Float64,4}
    m, k = size(A)
    kb, n = size(B)
    k == kb || throw(DimensionMismatch(
        "inner dimensions differ: size(A,2)=$k and size(B,1)=$kb"))
    size(C) == (m,n) || throw(DimensionMismatch(
        "C has size $(size(C)); expected ($m, $n)"))
    (Base.mightalias(C, A) || Base.mightalias(C, B)) && throw(ArgumentError(
        "gemm! does not permit C to alias A or B"))

    a = _coerce(M, α)
    if isone(a) && m >= 8 && n >= 2 && k >= 4
        return _gemm_vec2col!(C, A, B, Val(8); α=a, β=β)
    end
    return _gemm_streaming_kernel!(C, A, B; α=a, β=β)
end

end # @eval MFLinearAlgebra
