@testset "M6 Float64x5 reciprocal seed/correction A/B" begin
    Random.seed!(0x6ec5_2026)

    T = MultiFloat{Float64,5}
    recip5 = ExperimentalArithmetic.recip5_safe
    ref_inv = ExperimentalArithmetic.reference_inv

    function pow2_recip5(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    function dense_recip5(; words=5, emin=-120, emax=120)
        numerator = BigInt(0)
        for _ in 1:words
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        numerator |= BigInt(1) << (64words - 1)
        q = numerator // (BigInt(1) << 64words)
        q *= pow2_recip5(rand(emin:emax))
        q = rand(Bool) ? q : -q
        return setprecision(BigFloat, 1536) do
            T(BigFloat(q))
        end
    end

    exact_residual(x::T, z::T) =
        one(Rational{BigInt}) - Rational{BigInt}(x) * Rational{BigInt}(z)

    function check_case(x::T)
        oracle = ref_inv(x)
        x4seed = ExperimentalArithmetic._recip5_seed_x4(x)
        x4one = ExperimentalArithmetic._recip5_x4_one_correction(x)
        f64three = ExperimentalArithmetic._recip5_float64_three_corrections(x)
        z = recip5(x)

        @test MultiFloats.isnormalized(x4seed)
        @test MultiFloats.isnormalized(x4one)
        @test MultiFloats.isnormalized(f64three)
        @test z === x4one
        @test z === oracle
        @test f64three === oracle

        r0 = abs(exact_residual(x, x4seed))
        r1 = abs(exact_residual(x, x4one))
        @test r1 <= r0
        setprecision(BigFloat, 1024) do
            b0 = BigFloat(r0)
            b1 = BigFloat(r1)
            @test b1 <= max(BigFloat(8) * b0^2, BigFloat(2)^(-260))
        end

        qref = inv(Rational{BigInt}(x))
        err = abs(Rational{BigInt}(z) - qref)
        setprecision(BigFloat, 1024) do
            @test BigFloat(err) / abs(BigFloat(qref)) <= BigFloat(2)^(-250)
        end
        return z
    end

    @testset "ordinary and scaled inputs" begin
        for _ in 1:100
            check_case(dense_recip5(emin=-40, emax=40))
        end
        for _ in 1:100
            check_case(dense_recip5(emin=-120, emax=120))
        end
    end

    @testset "identities, signs, and powers of two" begin
        o = one(T)
        @test recip5(o) === o
        @test recip5(T(-1.0)) === T(-1.0)
        for e in (-200, -100, 0, 100, 200)
            x = T(BigFloat(pow2_recip5(e)))
            nx = T(BigFloat(-pow2_recip5(e)))
            @test recip5(x) === T(BigFloat(pow2_recip5(-e)))
            @test recip5(nx) === T(BigFloat(-pow2_recip5(-e)))
        end
    end

    @testset "explicit candidate domain" begin
        @test_throws DomainError recip5(T(0.0))
        @test_throws DomainError recip5(T(Inf))

        unnormalized = T((1.0, 1.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError recip5(unnormalized)

        tiny_normal = T(floatmin(Float64))
        @test recip5(tiny_normal) === ref_inv(tiny_normal)

        tiny_subnormal = T(nextfloat(0.0))
        @test_throws DomainError recip5(tiny_subnormal)

        huge = T(floatmax(Float64))
        @test_throws DomainError recip5(huge)
    end
end
