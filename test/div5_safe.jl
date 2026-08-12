@testset "M6 Float64x5 direct-residual division candidate" begin
    Random.seed!(0x6d15_2026)

    T = MultiFloat{Float64,5}
    div5 = ExperimentalArithmetic.div5_safe
    ref_div = ExperimentalArithmetic.reference_div

    function pow2_div5(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    function dense_div5(; words=5, emin=-80, emax=80)
        numerator = BigInt(0)
        for _ in 1:words
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        numerator |= BigInt(1) << (64words - 1)
        q = numerator // (BigInt(1) << 64words)
        q *= pow2_div5(rand(emin:emax))
        q = rand(Bool) ? q : -q
        return setprecision(BigFloat, 1536) do
            T(BigFloat(q))
        end
    end

    function check_case(x::T, y::T)
        q0, invy = ExperimentalArithmetic._div5_seed(x, y)
        q1, residual = ExperimentalArithmetic._div5_correct_once(x, y, q0, invy)
        q = div5(x, y)
        oracle = ref_div(x, y)

        @test MultiFloats.isnormalized(q0)
        @test MultiFloats.isnormalized(invy)
        @test MultiFloats.isnormalized(residual)
        @test MultiFloats.isnormalized(q1)
        @test q === q1

        # Sole numerical acceptance gate: the corrected quotient must equal the
        # independent adaptive division oracle bit-for-bit.
        @test q === oracle

        # The reported residual is itself one accepted direct FMA result.
        nq0 = ExperimentalArithmetic._neg_recip5_safe(q0)
        @test residual === ExperimentalArithmetic.fma5_safe(nq0, y, x)
        return nothing
    end

    @testset "ordinary and scaled inputs" begin
        for _ in 1:80
            check_case(dense_div5(emin=-40, emax=40), dense_div5(emin=-40, emax=40))
        end
        for _ in 1:80
            check_case(dense_div5(emin=-100, emax=100), dense_div5(emin=-100, emax=100))
        end
    end

    @testset "identities and signs" begin
        o = one(T)
        z = zero(T)
        for _ in 1:30
            x = dense_div5(emin=-80, emax=80)
            @test div5(x, o) === x
            @test div5(x, x) === o
            @test div5(z, x) === z

            nx = setprecision(BigFloat, 1536) do
                T(BigFloat(-Rational{BigInt}(x)))
            end
            @test div5(nx, x) == -o
            @test div5(x, nx) == -o
        end
    end

    @testset "powers of two" begin
        for ex in (-120, -60, 0, 60, 120), ey in (-120, -60, 0, 60, 120)
            x = T(BigFloat(pow2_div5(ex)))
            y = T(BigFloat(pow2_div5(ey)))
            @test div5(x, y) === T(BigFloat(pow2_div5(ex - ey)))
        end
    end

    @testset "explicit domain" begin
        @test_throws DomainError div5(T(1.0), T(0.0))
        @test_throws DomainError div5(T(Inf), T(1.0))
        @test_throws DomainError div5(T(1.0), T(Inf))

        unnormalized = T((1.0, 1.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError div5(unnormalized, T(1.0))
        @test_throws ArgumentError div5(T(1.0), unnormalized)

        tiny_subnormal = T(nextfloat(0.0))
        @test_throws DomainError div5(T(1.0), tiny_subnormal)

        huge = T(floatmax(Float64))
        @test_throws DomainError div5(T(1.0), huge)
    end
end
