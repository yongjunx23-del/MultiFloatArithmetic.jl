@testset "safe Float64x5 direct FMA candidate" begin
    Random.seed!(0xf5a5_2026)

    T = MultiFloat{Float64,5}
    fma5 = ExperimentalArithmetic.fma5_safe
    reference_fma = ExperimentalArithmetic.reference_fma

    function pow2_fma5(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    function dense_fma5(; words=5, emin=-70, emax=70)
        numerator = BigInt(0)
        for _ in 1:words
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        numerator |= BigInt(1) << (64words - 1)
        q = numerator // (BigInt(1) << 64words)
        q *= pow2_fma5(rand(emin:emax))
        q = rand(Bool) ? q : -q
        return setprecision(BigFloat, 1536) do
            T(BigFloat(q))
        end
    end

    function exact_fma5_tail(full)
        total = zero(Rational{BigInt})
        for i in 6:55
            total += Rational{BigInt}(full[i])
        end
        return total
    end

    function check_fma5_case(x::T, y::T, c::T)
        for i in 1:5, j in 1:5
            xi = x._limbs[i]
            yj = y._limbs[j]
            p, e = MultiFloats.two_prod(xi, yj)
            @test Rational{BigInt}(p) + Rational{BigInt}(e) ==
                  Rational{BigInt}(xi) * Rational{BigInt}(yj)
        end

        z = fma5(x, y, c)
        terms = ExperimentalArithmetic._fma5_terms(x._limbs, y._limbs, c._limbs)
        terms_swapped = ExperimentalArithmetic._fma5_terms(y._limbs, x._limbs, c._limbs)
        full = ExperimentalArithmetic._fma5_full_limbs(x._limbs, y._limbs, c._limbs)

        exact = Rational{BigInt}(x) * Rational{BigInt}(y) + Rational{BigInt}(c)
        discarded = exact_fma5_tail(full)

        @test length(full) == 55
        @test terms === terms_swapped
        @test MultiFloats.isnormalized(full)
        @test MultiFloats.isnormalized(z)
        @test Rational{BigInt}(z) + discarded == exact
        @test z === reference_fma(x, y, c)
        @test z === fma5(y, x, c)
        return z
    end

    @testset "dense ordinary/scaled" begin
        for _ in 1:100
            check_fma5_case(dense_fma5(), dense_fma5(), dense_fma5())
        end
        for _ in 1:100
            check_fma5_case(
                dense_fma5(emin=-110, emax=110),
                dense_fma5(emin=-110, emax=110),
                dense_fma5(emin=-110, emax=110),
            )
        end
    end

    @testset "exact identities" begin
        z = zero(T)
        o = one(T)
        for _ in 1:30
            x = dense_fma5(emin=-80, emax=80)
            c = dense_fma5(emin=-80, emax=80)
            @test fma5(x, o, z) === x
            @test fma5(x, z, c) === c
            @test fma5(z, x, c) === c

            nx = setprecision(BigFloat, 1536) do
                T(BigFloat(-Rational{BigInt}(x)))
            end
            @test iszero(fma5(x, o, nx))
            check_fma5_case(x, o, nx)
        end
    end

    @testset "destructive cancellation" begin
        for depth in (80, 160, 230, 270, 320)
            for _ in 1:30
                x = dense_fma5(emin=-50, emax=50)
                y = dense_fma5(emin=-50, emax=50)
                qxy = Rational{BigInt}(x) * Rational{BigInt}(y)
                lead_exp = exponent(x._limbs[1]) + exponent(y._limbs[1])
                delta = pow2_fma5(lead_exp - depth)
                c = setprecision(BigFloat, 1536) do
                    T(BigFloat(-qxy + (rand(Bool) ? delta : -delta)))
                end
                check_fma5_case(x, y, c)
            end
        end
    end

    @testset "power-of-two and boundary cases" begin
        for ex in (-100, 0, 100), ey in (-100, 0, 100), ec in (-100, 0, 100)
            x = T(BigFloat(pow2_fma5(ex)))
            y = T(BigFloat(pow2_fma5(ey)))
            c = T(BigFloat(pow2_fma5(ec)))
            check_fma5_case(x, y, c)
            check_fma5_case(x, y, T(BigFloat(-pow2_fma5(ec))))
        end
    end

    @testset "explicit domain" begin
        @test_throws DomainError fma5(T(Inf), T(1.0), T(0.0))
        @test_throws DomainError fma5(T(1.0), T(1.0), T(Inf))
        unnormalized = T((1.0, 1.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError fma5(unnormalized, T(1.0), T(0.0))

        huge = T(ldexp(1.0, 800))
        large = T(ldexp(1.0, 400))
        @test_throws DomainError fma5(huge, large, T(0.0))

        tiny = T(ldexp(1.0, -700))
        small = T(ldexp(1.0, -400))
        @test_throws DomainError fma5(tiny, small, T(0.0))
    end
end

include("fma6_fma8_safe.jl")
include("recip5_safe.jl")
