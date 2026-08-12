using LinearAlgebra

const MFLA = MultiFloatArithmetic.MFLinearAlgebra

# Load the experimental scheduling layer only on this branch. Promotion will
# move the selected kernel into the package include graph after the A/B gate.
include(joinpath(@__DIR__, "..", "src", "linear_algebra_blocked.jl"))

function bigvec(x)
    setprecision(BigFloat, 1024) do
        [BigFloat(v) for v in x]
    end
end

function bigmat(A)
    setprecision(BigFloat, 1024) do
        [BigFloat(A[i,j]) for i in axes(A,1), j in axes(A,2)]
    end
end

function max_relative_error(actual, reference)
    setprecision(BigFloat, 1024) do
        worst = BigFloat(0)
        for idx in eachindex(actual, reference)
            a = BigFloat(actual[idx])
            r = BigFloat(reference[idx])
            scale = max(abs(r), one(BigFloat))
            worst = max(worst, abs(a-r) / scale)
        end
        return worst
    end
end

function la_tolerance(::Type{MultiFloat{T,N}}) where {T,N}
    # Linear-algebra tests exercise repeated FMA chains rather than one kernel.
    # Leave ~28 bits of headroom for accumulation length, scaling, and the
    # independently rounded BigFloat -> MultiFloat test inputs.
    return setprecision(BigFloat, 1024) do
        BigFloat(2)^(-(precision(T)*N - 28))
    end
end

function random_mf_matrix(::Type{M}, m, n; scale=1.0) where {M<:MultiFloat}
    [M(scale * (2rand() - 1)) for _ in 1:m, _ in 1:n]
end

function random_mf_vector(::Type{M}, n; scale=1.0) where {M<:MultiFloat}
    [M(scale * (2rand() - 1)) for _ in 1:n]
end

@testset "MFLinearAlgebra API and FMA chaining" begin
    Random.seed!(0x1a_2026)

    # Cover every public hot width/type combination claimed by MFLinearAlgebra.
    for M in (Float32x2, Float32x3, Float32x4,
              Float64x2, Float64x3, Float64x4)
        tol = la_tolerance(M)

        @testset "$(M) dot/axpy" begin
            x = random_mf_vector(M, 17)
            y = random_mf_vector(M, 17)
            xb, yb = bigvec(x), bigvec(y)

            d = MFLA.mfdot(x, y)
            dref = setprecision(BigFloat, 1024) do
                sum(xb[i] * yb[i] for i in eachindex(xb))
            end
            @test abs(BigFloat(d) - dref) <= tol * max(abs(dref), one(BigFloat))
            @test MultiFloats.isnormalized(d)

            y0 = copy(y)
            α = M(0.375)
            MFLA.axpy!(α, x, y0)
            yref = setprecision(BigFloat, 1024) do
                BigFloat(α) .* xb .+ yb
            end
            @test max_relative_error(y0, yref) <= 4tol
            @test all(MultiFloats.isnormalized, y0)
        end

        @testset "$(M) GEMV" begin
            A = random_mf_matrix(M, 9, 7)
            x = random_mf_vector(M, 7)
            y = random_mf_vector(M, 9)
            Ab, xb, yb = bigmat(A), bigvec(x), bigvec(y)

            yfast = MFLA.gemv(A, x)
            yref = setprecision(BigFloat, 1024) do
                Ab * xb
            end
            @test max_relative_error(yfast, yref) <= 8tol
            @test all(MultiFloats.isnormalized, yfast)

            α = M(0.75)
            β = M(-0.25)
            yscaled = copy(y)
            MFLA.gemv!(yscaled, A, x; α=α, β=β)
            yscaled_ref = setprecision(BigFloat, 1024) do
                BigFloat(α) .* (Ab * xb) .+ BigFloat(β) .* yb
            end
            @test max_relative_error(yscaled, yscaled_ref) <= 12tol
            @test all(MultiFloats.isnormalized, yscaled)
        end

        @testset "$(M) GEMM" begin
            A = random_mf_matrix(M, 8, 6)
            B = random_mf_matrix(M, 6, 7)
            Ab, Bb = bigmat(A), bigmat(B)

            C = MFLA.gemm(A, B)
            Cref = setprecision(BigFloat, 1024) do
                Ab * Bb
            end
            @test max_relative_error(C, Cref) <= 12tol
            @test all(MultiFloats.isnormalized, C)

            C0 = random_mf_matrix(M, 8, 7)
            C0b = bigmat(C0)
            α = M(0.5)
            β = M(0.125)
            Cscaled = copy(C0)
            MFLA.gemm!(Cscaled, A, B; α=α, β=β)
            Cscaled_ref = setprecision(BigFloat, 1024) do
                BigFloat(α) .* (Ab * Bb) .+ BigFloat(β) .* C0b
            end
            @test max_relative_error(Cscaled, Cscaled_ref) <= 16tol
            @test all(MultiFloats.isnormalized, Cscaled)

            @test_throws DimensionMismatch MFLA.gemm!(zeros(M, 2, 2), A, B)
            @test_throws ArgumentError MFLA.gemm!(A, A, Matrix{M}(I, 6, 6))
        end

        @testset "$(M) SYRK" begin
            A = random_mf_matrix(M, 8, 5)
            Ab = bigmat(A)
            C = MFLA.syrk(A; uplo=:L)
            Cref = setprecision(BigFloat, 1024) do
                Ab * transpose(Ab)
            end
            @test max_relative_error(C, Cref) <= 14tol
            @test C == transpose(C)
            @test all(MultiFloats.isnormalized, C)
        end

        @testset "$(M) triangular solve" begin
            n = 7
            L = zeros(M, n, n)
            @inbounds for j in 1:n
                L[j,j] = M(2.0 + 0.1j)
                for i in j+1:n
                    L[i,j] = M(0.05 * (2rand() - 1))
                end
            end
            Lb = bigmat(L)

            xtrue = random_mf_vector(M, n)
            bbig = setprecision(BigFloat, 1024) do
                Lb * bigvec(xtrue)
            end
            b = M.(bbig)
            MFLA.trsv!(L, b; uplo=:L)
            @test max_relative_error(b, bigvec(xtrue)) <= 20tol
            @test all(MultiFloats.isnormalized, b)

            Xtrue = random_mf_matrix(M, n, 3)
            Bbig = setprecision(BigFloat, 1024) do
                Lb * bigmat(Xtrue)
            end
            B = M.(Bbig)
            MFLA.trsm!(L, B; uplo=:L)
            @test max_relative_error(B, bigmat(Xtrue)) <= 24tol
            @test all(MultiFloats.isnormalized, B)
        end

        @testset "$(M) Cholesky" begin
            n = 6
            Rbig = setprecision(BigFloat, 1024) do
                [BigFloat(0.25 * (2rand() - 1)) for _ in 1:n, _ in 1:n]
            end
            Abig = setprecision(BigFloat, 1024) do
                Rbig * transpose(Rbig) + BigFloat(2) * Matrix{BigFloat}(I, n, n)
            end
            A = M.(Abig)
            L = copy(A)
            MFLA.potrf!(L; uplo=:L)

            @test all(iszero, triu(L, 1))
            @test all(MultiFloats.isnormalized, L)
            Lb = bigmat(L)
            reconstructed = setprecision(BigFloat, 1024) do
                Lb * transpose(Lb)
            end
            @test max_relative_error(reconstructed, Abig) <= 32tol
        end
    end
end

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

@testset "Float64x4 blocked Cholesky preserves the scalar schedule" begin
    for n in (17, 33)
        A = blocked_test_spd_x4(n)
        baseline = copy(A)
        MFLA._potrf_unblocked_kernel!(baseline)
        @test all(MultiFloats.isnormalized, baseline)

        for bs in (8, 16)
            blocked = copy(A)
            MFLA._potrf_blocked_vec8!(blocked, Val(bs))
            @test all(i -> blocked[i] === baseline[i], eachindex(blocked, baseline))
            @test all(MultiFloats.isnormalized, blocked)
        end
    end
end
