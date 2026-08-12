# Explicit MultiFloatVec microkernels layered on top of the auditable scalar
# MFLinearAlgebra baseline.  Keep these implementation details in the existing
# MFLinearAlgebra namespace so they can be A/B tested without widening the
# public API.

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

"""
    _gemm_vec_accumulate!(C, A, B, Val(W))

Internal explicit-SIMD GEMM accumulation experiment for `α = 1` after `C` has
already been initialized/scaled for `β`.

The kernel keeps W independent output rows in one `MultiFloatVec` accumulator
across the complete k reduction.  Each lane sees the same `p = 1:k` reduction
order as the scalar streaming kernel, so a successful implementation must be
bitwise identical lane by lane.  A scalar remainder handles rows not divisible
by W.
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
                bv = V(B[p,j])
                cv = fma_fast(av, bv, cv)
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
    _gemm_vec!(C, A, B, Val(W); α=1, β=0)

Full-semantics wrapper for the explicit-SIMD experiment.  The current
microkernel intentionally specializes only the dominant `α == 1` route; other
scaling cases fall back to the public scalar-streaming implementation so the
benchmark cannot manufacture a speedup by changing semantics.
"""
function _gemm_vec!(
    C::AbstractMatrix{MultiFloat{T,N}},
    A::AbstractMatrix{MultiFloat{T,N}},
    B::AbstractMatrix{MultiFloat{T,N}},
    width::Val{W};
    α=one(MultiFloat{T,N}),
    β=zero(MultiFloat{T,N}),
) where {T,N,W}
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
    if !isone(a)
        return gemm!(C, A, B; α=a, β=b)
    end

    if iszero(b)
        fill!(C, zero(M))
    elseif !isone(b)
        @inbounds for j in 1:n
            @simd for i in 1:m
                C[i,j] = _scale(C[i,j], b)
            end
        end
    end

    return _gemm_vec_accumulate!(C, A, B, width)
end

end # @eval MFLinearAlgebra
