using Random
using Test

const DenseGEMMM = Float64x4

function _dense_gemm_case(m, k, n; α=one(DenseGEMMM), β=zero(DenseGEMMM))
    A = rand(DenseGEMMM, m, k)
    B = rand(DenseGEMMM, k, n)
    C0 = rand(DenseGEMMM, m, n)
    Cfast = copy(C0)
    Cstream = copy(C0)

    MFLA.gemm!(Cfast, A, B; α=α, β=β)
    MFLA._gemm_streaming_kernel!(Cstream, A, B; α=α, β=β)

    @test all(i -> Cfast[i] === Cstream[i], eachindex(Cfast))
    @test all(MultiFloats.isnormalized, Cfast)
    return nothing
end

@testset "Float64x4 dense MR8xNR2 GEMM fast path" begin
    Random.seed!(0xd8e4_2026)

    # Fast-path shapes, including simultaneous row and column remainders.
    _dense_gemm_case(8, 4, 2)
    _dense_gemm_case(16, 9, 6)
    _dense_gemm_case(19, 11, 7)
    _dense_gemm_case(25, 13, 5; β=DenseGEMMM(0.25))
    _dense_gemm_case(24, 12, 8; β=one(DenseGEMMM))

    # Non-unit alpha deliberately remains on the streaming fallback.
    _dense_gemm_case(19, 11, 7; α=DenseGEMMM(0.5), β=DenseGEMMM(0.125))
    _dense_gemm_case(19, 11, 7; α=zero(DenseGEMMM), β=DenseGEMMM(-0.25))

    # Threshold boundaries must retain the same public semantics.
    _dense_gemm_case(7, 11, 7)
    _dense_gemm_case(19, 3, 7)
    _dense_gemm_case(19, 11, 1)
end
