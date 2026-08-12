# Blocked dense factorization kernels layered on top of MFLinearAlgebra.
# Each matrix element sees the same p=1,2,... chained-FMA order and operand
# roles as the unblocked baseline. Blocking only reschedules independent element
# updates; Vec8 only executes independent rows in parallel.

@eval MFLinearAlgebra begin

function _potrf_unblocked_kernel!(
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

function _chol_trailing_update_vec8!(
    A::Matrix{MultiFloat{Float64,4}},
    row0::Int,
    p0::Int,
    p1::Int,
)
    V = MultiFloatVec{8,Float64,4}
    n = size(A,1)
    row0 > n && return A

    @inbounds for p in p0:p1
        for j in row0:n
            coeff = A[j,p]
            nrows = n - j + 1
            full = nrows - rem(nrows, 8)
            stop = j + full - 1
            for i0 in j:8:stop
                cv = _load_mfvec(V, A, i0, j)
                av = _load_mfvec(V, A, i0, p)
                cv = fma_fast(-av, V(coeff), cv)
                _store_mfvec!(A, i0, j, cv)
            end
            for i in stop + 1:n
                A[i,j] = fma_fast(-A[i,p], coeff, A[i,j])
            end
        end
    end
    return A
end

function _chol_factor_panel_vec8!(
    A::Matrix{MultiFloat{Float64,4}},
    j0::Int,
    j1::Int,
)
    V = MultiFloatVec{8,Float64,4}
    n = size(A,1)

    @inbounds for j in j0:j1
        d = A[j,j]
        for p in j0:j-1
            d = fma_fast(-A[j,p], A[j,p], d)
        end
        A[j,j] = sqrt_r(d)
        diag = A[j,j]

        firstrow = j + 1
        firstrow > n && continue
        nrows = n - firstrow + 1
        full = nrows - rem(nrows, 8)
        stop = firstrow + full - 1

        for i0 in firstrow:8:stop
            sv = _load_mfvec(V, A, i0, j)
            for p in j0:j-1
                av = _load_mfvec(V, A, i0, p)
                sv = fma_fast(-av, V(A[j,p]), sv)
            end
            for lane in 1:8
                A[i0 + lane - 1,j] = div_r(sv[lane], diag)
            end
        end
        for i in stop + 1:n
            s = A[i,j]
            for p in j0:j-1
                s = fma_fast(-A[i,p], A[j,p], s)
            end
            A[i,j] = div_r(s, diag)
        end
    end
    return A
end

"""
    _potrf_blocked_vec8!(A, Val(BS); clean=true)

Lower Cholesky kernel for dense Float64x4 matrices. It is deliberately
bitwise-equivalent to `_potrf_unblocked_kernel!` and uses explicit Vec8 only
across independent rows.
"""
function _potrf_blocked_vec8!(
    A::Matrix{MultiFloat{Float64,4}},
    ::Val{BS};
    clean::Bool=true,
) where {BS}
    BS > 0 || throw(ArgumentError("Cholesky block size must be positive"))
    n, n2 = size(A)
    n == n2 || throw(DimensionMismatch("Cholesky factorization requires a square matrix"))

    j0 = 1
    while j0 <= n
        j1 = min(n, j0 + BS - 1)
        _chol_factor_panel_vec8!(A, j0, j1)
        row0 = j1 + 1
        if row0 <= n
            _chol_trailing_update_vec8!(A, row0, j0, j1)
        end
        j0 = row0
    end

    if clean
        M = MultiFloat{Float64,4}
        @inbounds for j in 1:n, i in 1:j-1
            A[i,j] = zero(M)
        end
    end
    return A
end

# Hosted Zen3 A/B, after whole-matrix bitwise-equality gates:
#   n=32: 1.70x, n=64: 2.34x, n=96: 2.74x (BS=8 versus unblocked).
# Keep the dispatch narrower than the implementation: only dense Float64x4,
# lower factorization, and the smallest size actually covered by the timing A/B.
function potrf!(
    A::Matrix{MultiFloat{Float64,4}};
    uplo::Symbol=:L,
    clean::Bool=true,
)
    n = size(A,1)
    if uplo === :L && n >= 32
        return _potrf_blocked_vec8!(A, Val(8); clean=clean)
    end
    return _potrf_unblocked_kernel!(A; uplo=uplo, clean=clean)
end

end # @eval MFLinearAlgebra
