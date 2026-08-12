@testset "M6 Float64x8 direct-FMA reciprocal A/B" begin
    Random.seed!(0x6ec8_2026)

    T = MultiFloat{Float64,8}
    inv8 = ExperimentalArithmetic.inv8_safe
    ref_inv = ExperimentalArithmetic.reference_inv

    function pow2_inv8(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    function dense_inv8(; words=8, emin=-80, emax=80)
        numerator = BigInt(0)
        for _ in 1:words
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        numerator |= BigInt(1) << (64words - 1)
        q = numerator // (BigInt(1) << 64words)
        q *= pow2_inv8(rand(emin:emax))
        q = rand(Bool) ? q : -q
        return setprecision(BigFloat, 2048) do
            T(BigFloat(q))
        end
    end

    exact_residual(x::T, y::T) =
        one(Rational{BigInt}) - Rational{BigInt}(x) * Rational{BigInt}(y)

    function check_case(x::T)
        oracle = ref_inv(x)
        seed = ExperimentalArithmetic._inv8_seed_x4(x)
        one_step = ExperimentalArithmetic._inv8_direct_one_correction(x)
        two_step = ExperimentalArithmetic._inv8_direct_two_corrections(x)
        three_step = ExperimentalArithmetic._inv8_direct_three_corrections(x)
        public = inv8(x)

        @test MultiFloats.isnormalized(seed)
        @test MultiFloats.isnormalized(one_step)
        @test MultiFloats.isnormalized(two_step)
        @test MultiFloats.isnormalized(three_step)
        @test public === three_step

        # Sole numerical acceptance gate. One/two-step equality is diagnostic;
        # three steps must equal the independent adaptive oracle bit-for-bit.
        @test three_step === oracle

        r0 = abs(exact_residual(x, seed))
        r1 = abs(exact_residual(x, one_step))
        r2 = abs(exact_residual(x, two_step))
        r3 = abs(exact_residual(x, three_step))
        @test r1 <= r0
        @test r2 <= r1
        @test r3 <= r2
        return nothing
    end

    @testset "ordinary and scaled inputs" begin
        for _ in 1:20
            check_case(dense_inv8(emin=-40, emax=40))
        end
        for _ in 1:20
            check_case(dense_inv8(emin=-100, emax=100))
        end
    end

    @testset "identities, signs, and powers of two" begin
        o = one(T)
        @test inv8(o) === o
        @test inv8(T(-1.0)) === T(-1.0)
        for e in (-160, -80, 0, 80, 160)
            x = T(BigFloat(pow2_inv8(e)))
            nx = T(BigFloat(-pow2_inv8(e)))
            @test inv8(x) === T(BigFloat(pow2_inv8(-e)))
            @test inv8(nx) === T(BigFloat(-pow2_inv8(-e)))
        end
    end

    @testset "explicit candidate domain" begin
        @test_throws DomainError inv8(T(0.0))
        @test_throws DomainError inv8(T(Inf))

        unnormalized = T((1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError inv8(unnormalized)

        tiny_normal = T(floatmin(Float64))
        @test inv8(tiny_normal) === ref_inv(tiny_normal)

        tiny_subnormal = T(nextfloat(0.0))
        @test_throws DomainError inv8(tiny_subnormal)

        huge = T(floatmax(Float64))
        @test_throws DomainError inv8(huge)
    end
end
