@testset "safe Float64x5 addition candidate" begin
    Random.seed!(0xad05_2026)

    T = MultiFloat{Float64,5}
    add5 = ExperimentalArithmetic.add5_safe
    reference_add = ExperimentalArithmetic.reference_add

    function exact_limb_sum(limbs)
        total = zero(Rational{BigInt})
        for limb in limbs
            total += Rational{BigInt}(limb)
        end
        return total
    end

    function pow2_rational(e::Int)
        if e >= 0
            return (BigInt(1) << e) // BigInt(1)
        end
        return BigInt(1) // (BigInt(1) << (-e))
    end

    function check_add5_case(x::T, y::T)
        z = add5(x, y)
        full = ExperimentalArithmetic._add_safe_full_limbs(x._limbs, y._limbs)
        tail = ntuple(i -> full[i + 5], Val{5}())

        exact = Rational{BigInt}(x) + Rational{BigInt}(y)
        represented = Rational{BigInt}(z)
        discarded = exact_limb_sum(tail)

        @test MultiFloats.isnormalized(z)
        @test MultiFloats.isnormalized(full)
        @test represented + discarded == exact
        @test z === reference_add(x, y)
        @test z === add5(y, x)
        return z
    end

    @testset "ordinary and wide-exponent inputs" begin
        for _ in 1:250
            x = rand(T)
            y = rand(T)
            check_add5_case(x, y)
        end

        for _ in 1:250
            x = wide_rand(T; emin=-400, emax=400)
            y = wide_rand(T; emin=-400, emax=400)
            check_add5_case(x, y)
        end
    end

    @testset "exact identities" begin
        for _ in 1:50
            x = wide_rand(T; emin=-300, emax=300)
            @test add5(x, zero(T)) === x
            @test add5(zero(T), x) === x
            @test iszero(add5(x, -x))
            check_add5_case(x, -x)
        end
    end

    @testset "near cancellation" begin
        for bits in (100, 180, 240, 280)
            for _ in 1:40
                x = wide_rand(T; emin=-200, emax=200)
                qx = Rational{BigInt}(x)
                lead = x._limbs[1]
                e = iszero(lead) ? 0 : exponent(lead)
                delta = pow2_rational(e - bits)
                y = T(-qx + (rand(Bool) ? delta : -delta))
                check_add5_case(x, y)
            end
        end
    end

    @testset "power-of-two and carry boundaries" begin
        for e in (-300, -50, -1, 0, 1, 50, 300)
            base = pow2_rational(e)
            for depth in (53, 106, 159, 212, 265, 270, 300)
                delta = pow2_rational(e - depth)

                x = T(base - delta)
                y = T(delta)
                z = check_add5_case(x, y)
                if Rational{BigInt}(x) + Rational{BigInt}(y) == base
                    @test z === T(base)
                end

                x = T(base + delta)
                y = T(-delta)
                z = check_add5_case(x, y)
                if Rational{BigInt}(x) + Rational{BigInt}(y) == base
                    @test z === T(base)
                end

                check_add5_case(T(base), T(delta))
                check_add5_case(T(-base), T(delta))
            end
        end
    end

    @testset "explicit input contract" begin
        @test_throws DomainError add5(T(Inf), T(1.0))
        unnormalized = T((1.0, 1.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError add5(unnormalized, T(1.0))
    end
end
