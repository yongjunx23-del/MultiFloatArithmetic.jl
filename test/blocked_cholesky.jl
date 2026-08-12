const BlockedLA = MultiFloatArithmetic.MFLinearAlgebra

function blocked_test_spd_x4(n)
    M = Float64x4
    Random.seed!(0xbc40 + n)
    L = zeros(M, n, n)
    @inbounds for j in 1:n
        L[j,j] = M(1.4 + 0.01j)
        for i in j+1:n
            L[i,j] = M(0.02 * (2rand() - 1))
        end
    end

    A = zeros(M, n, n)
    @inbounds for j in 1:n
        for i in j:n
            s = zero(M)
            for p in 1:j
                s = fma_fast(L[i,p], L[j,p], s)
            end
            A[i,j] = s
            A[j,i] = s
        end
    end
    return A
end

@testset "Float64x4 blocked Cholesky preserves scalar schedule" begin
    for n in (17, 33)
        A = blocked_test_spd_x4(n)
        baseline = copy(A)
        BlockedLA._potrf_unblocked_kernel!(baseline)
        @test all(MultiFloats.isnormalized, baseline)

        public = copy(A)
        BlockedLA.potrf!(public; uplo=:L)
        @test all(i -> public[i] === baseline[i], eachindex(public, baseline))
        @test all(MultiFloats.isnormalized, public)

        for bs in (8, 16)
            blocked = copy(A)
            BlockedLA._potrf_blocked_vec8!(blocked, Val(bs))
            @test all(i -> blocked[i] === baseline[i], eachindex(blocked, baseline))
            @test all(MultiFloats.isnormalized, blocked)
        end
    end

    # The optimized dispatch is lower-only. A large upper factorization must
    # still take the unblocked fallback without changing semantics.
    A = blocked_test_spd_x4(33)
    upper_ref = copy(A)
    upper_public = copy(A)
    BlockedLA._potrf_unblocked_kernel!(upper_ref; uplo=:U)
    BlockedLA.potrf!(upper_public; uplo=:U)
    @test all(i -> upper_public[i] === upper_ref[i], eachindex(upper_ref, upper_public))
    @test all(MultiFloats.isnormalized, upper_public)
end
