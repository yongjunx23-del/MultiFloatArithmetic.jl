@testset "safe Float64x5 multiplication candidate" begin
    Random.seed!(0x4d55_2026)

    T = MultiFloat{Float64,5}
    mul5 = ExperimentalArithmetic.mul5_safe
    reference_mul = ExperimentalArithmetic.reference_mul

    function pow2_mul5(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    function dense_mul5(; words=5, emin=-80, emax=80)
        numerator = BigInt(0)
        for _ in 1:words
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        numerator |= BigInt(1) << (64words - 1)
        q = numerator // (BigInt(1) << 64words)
        e = rand(emin:emax)
        q *= pow2_mul5(e)
        q = rand(Bool) ? q : -q
        return setprecision(BigFloat, 1024) do
            T(BigFloat(q))
        end
    end

    function exact_tail_sum(full)
        total = zero(Rational{BigInt})
        for i in 6:50
            total += Rational{BigInt}(full[i])
        end
        return total
    end

    function check_mul5_case(x::T, y::T)
        z = mul5(x, y)
        terms = ExperimentalArithmetic._mul5_terms(x._limbs, y._limbs)
        reverse_terms = ExperimentalArithmetic._mul5_terms(y._limbs, x._limbs)
        full = ExperimentalArithmetic._mul5_full_limbs(x._limbs, y._limbs)

        exact = Rational{BigInt}(x) * Rational{BigInt}(y)
        represented = Rational{BigInt}(z)
        discarded = exact_tail_sum(full)

        @test terms === reverse_terms
        @test MultiFloats.isnormalized(full)
        @test MultiFloats.isnormalized(z)
        @test represented + discarded == exact
        @test z === reference_mul(x, y)
        @test z === mul5(y, x)
        return z
    end

    @testset "dense ordinary and scaled inputs" begin
        for _ in 1:120
            check_mul5_case(dense_mul5(), dense_mul5())
        end

        for _ in 1:120
            x = dense_mul5(emin=-120, emax=120)
            y = dense_mul5(emin=-120, emax=120)
            check_mul5_case(x, y)
        end
    end

    @testset "exact identities and signs" begin
        z = zero(T)
        o = one(T)
        for _ in 1:40
            x = dense_mul5(emin=-100, emax=100)
            @test mul5(x, o) === x
            @test mul5(o, x) === x
            @test iszero(mul5(x, z))
            @test iszero(mul5(z, x))

            nx = setprecision(BigFloat, 1024) do
                T(BigFloat(-Rational{BigInt}(x)))
            end
            @test MultiFloats.isnormalized(nx)
            @test mul5(nx, o) === nx
            check_mul5_case(x, o)
            check_mul5_case(x, nx)
        end
    end

    @testset "power-of-two and near-boundary products" begin
        for ex in (-150, -50, 0, 50, 150)
            for ey in (-150, -50, 0, 50, 150)
                qx = pow2_mul5(ex)
                qy = pow2_mul5(ey)
                x = T(BigFloat(qx))
                y = T(BigFloat(qy))
                @test check_mul5_case(x, y) === T(BigFloat(qx * qy))
            end
        end

        for e in (-100, 0, 100)
            base = pow2_mul5(e)
            for depth in (53, 106, 159, 212, 250)
                delta = pow2_mul5(e - depth)
                x = T(BigFloat(base + delta))
                y = T(BigFloat(base - delta))
                check_mul5_case(x, y)
                check_mul5_case(x, T(BigFloat(-base + delta)))
            end
        end
    end

    @testset "explicit domain" begin
        @test_throws DomainError mul5(T(Inf), T(1.0))
        unnormalized = T((1.0, 1.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError mul5(unnormalized, T(1.0))

        # Leading pair product overflow.
        huge = T(ldexp(1.0, 800))
        large = T(ldexp(1.0, 400))
        @test_throws DomainError mul5(huge, large)

        # Pair product underflow to zero is outside the current M4 baseline.
        tiny = T(ldexp(1.0, -700))
        small = T(ldexp(1.0, -400))
        @test_throws DomainError mul5(tiny, small)
    end
end
