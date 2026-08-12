@testset "safe Float64x8 reciprocal candidate" begin
    Random.seed!(0x1a88_2026)

    T = MultiFloat{Float64,8}
    inv8 = ExperimentalArithmetic.inv8_safe
    ref_inv = ExperimentalArithmetic.reference_inv

    function repack_inv8(q::Rational{BigInt})
        return setprecision(BigFloat, 2048) do
            T(BigFloat(q))
        end
    end

    function moderate_inv8(; emin=-90, emax=90)
        q = Rational{BigInt}(rand(T))
        e = rand(emin:emax)
        if e >= 0
            q *= BigInt(1) << e
        else
            q /= BigInt(1) << (-e)
        end
        q = rand(Bool) ? q : -q
        x = repack_inv8(q)
        return iszero(x) ? T(1.25) : x
    end

    @testset "oracle agreement and normalization" begin
        for _ in 1:60
            x = moderate_inv8()
            z = inv8(x)
            @test MultiFloats.isnormalized(z)
            @test z === ref_inv(x)
        end
    end

    @testset "trusted seed and one-correction diagnostic" begin
        for _ in 1:30
            x = moderate_inv8(emin=-60, emax=60)
            seed = ExperimentalArithmetic._inv8_seed(x)
            one = ExperimentalArithmetic._inv8_safe_one_correction(x)
            @test MultiFloats.isnormalized(seed)
            @test MultiFloats.isnormalized(one)

            # The accepted two-correction baseline must agree with the oracle.
            @test inv8(x) === ref_inv(x)
        end
    end

    @testset "exact identities and powers of two" begin
        @test inv8(T(1.0)) === T(1.0)
        @test inv8(T(-1.0)) === T(-1.0)

        for e in (-200, -100, -1, 0, 1, 100, 200)
            x = T(ldexp(1.0, e))
            expected = T(ldexp(1.0, -e))
            @test inv8(x) === expected
            @test inv8(-T(ldexp(1.0, e))) === -expected
        end
    end

    @testset "safe negation helper" begin
        for _ in 1:20
            x = moderate_inv8(emin=-40, emax=40)
            nx = ExperimentalArithmetic._neg_safe(x)
            @test MultiFloats.isnormalized(nx)
            @test Rational{BigInt}(nx) == -Rational{BigInt}(x)
            @test ExperimentalArithmetic._neg_safe(nx) === x
        end
    end

    @testset "explicit domain" begin
        @test_throws DomainError inv8(zero(T))
        @test_throws DomainError inv8(T(Inf))
        unnormalized = T((1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError inv8(unnormalized)
    end
end
