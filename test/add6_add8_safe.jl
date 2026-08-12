@testset "safe Float64x6-x8 addition baselines" begin
    Random.seed!(0xad68_2026)

    reference_add = ExperimentalArithmetic.reference_add

    function exact_tail_sum(full, N)
        total = zero(Rational{BigInt})
        for i in N+1:2N
            total += Rational{BigInt}(full[i])
        end
        return total
    end

    function pow2q(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    for N in 6:8
        T = MultiFloat{Float64,N}
        addN = N == 6 ? ExperimentalArithmetic.add6_safe :
               N == 7 ? ExperimentalArithmetic.add7_safe :
                        ExperimentalArithmetic.add8_safe

        function check_case(x::T, y::T)
            z = addN(x, y)
            full = ExperimentalArithmetic._add_safe_full_limbs(x._limbs, y._limbs)
            exact = Rational{BigInt}(x) + Rational{BigInt}(y)
            tail = exact_tail_sum(full, N)

            @test MultiFloats.isnormalized(full)
            @test MultiFloats.isnormalized(z)
            @test Rational{BigInt}(z) + tail == exact
            @test z === reference_add(x, y)
            @test z === addN(y, x)
            @test z === ExperimentalArithmetic.add_safe(x, y)
            return z
        end

        @testset "Float64x$(N) ordinary/wide" begin
            for _ in 1:150
                check_case(rand(Bool) ? rand(T) : -rand(T),
                           rand(Bool) ? rand(T) : -rand(T))
            end
            for _ in 1:150
                check_case(wide_rand(T; emin=-350, emax=350),
                           wide_rand(T; emin=-350, emax=350))
            end
        end

        @testset "Float64x$(N) identities/cancellation" begin
            for _ in 1:25
                x = wide_rand(T; emin=-250, emax=250)
                @test addN(x, zero(T)) === x
                @test addN(zero(T), x) === x
                @test iszero(addN(x, -x))
                check_case(x, -x)
            end

            for depth in (100, 53N - 60, 53N - 10, 53N + 30)
                for _ in 1:20
                    x = wide_rand(T; emin=-150, emax=150)
                    qx = Rational{BigInt}(x)
                    lead = x._limbs[1]
                    e = iszero(lead) ? 0 : exponent(lead)
                    delta = pow2q(e - depth)
                    y = T(-qx + (rand(Bool) ? delta : -delta))
                    check_case(x, y)
                end
            end
        end

        @testset "Float64x$(N) power/carry boundaries" begin
            for e in (-200, 0, 200)
                base = pow2q(e)
                for depth in (53, 53(N - 1), 53N, 53N + 30)
                    delta = pow2q(e - depth)

                    x = T(base - delta)
                    y = T(delta)
                    z = check_case(x, y)
                    if Rational{BigInt}(x) + Rational{BigInt}(y) == base
                        @test z === T(base)
                    end

                    x = T(base + delta)
                    y = T(-delta)
                    z = check_case(x, y)
                    if Rational{BigInt}(x) + Rational{BigInt}(y) == base
                        @test z === T(base)
                    end

                    check_case(T(base), T(delta))
                    check_case(T(-base), T(delta))
                end
            end
        end
    end

    @testset "generic safe-add domain" begin
        T4 = MultiFloat{Float64,4}
        T9 = MultiFloat{Float64,9}
        @test_throws ArgumentError ExperimentalArithmetic.add_safe(T4(1.0), T4(2.0))
        @test_throws ArgumentError ExperimentalArithmetic.add_safe(T9(1.0), T9(2.0))
    end
end
