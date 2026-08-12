@testset "safe Float64x6-x8 multiplication baselines" begin
    Random.seed!(0x4d68_2026)

    reference_mul = ExperimentalArithmetic.reference_mul

    function pow2_mulN(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    function dense_mulN(::Type{T}, N; emin=-50, emax=50) where {T}
        numerator = BigInt(0)
        for _ in 1:N
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        numerator |= BigInt(1) << (64N - 1)
        q = numerator // (BigInt(1) << 64N)
        q *= pow2_mulN(rand(emin:emax))
        q = rand(Bool) ? q : -q
        return setprecision(BigFloat, 2048) do
            T(BigFloat(q))
        end
    end

    function exact_mul_tail(full, N)
        total = zero(Rational{BigInt})
        @inbounds for i in N+1:length(full)
            total += Rational{BigInt}(full[i])
        end
        return total
    end

    for N in 6:8
        T = MultiFloat{Float64,N}
        mulN = N == 6 ? ExperimentalArithmetic.mul6_safe :
               N == 7 ? ExperimentalArithmetic.mul7_safe :
                        ExperimentalArithmetic.mul8_safe

        check_case = function (x, y)
            # Primitive exactness is part of the accepted domain, not assumed.
            @inbounds for i in 1:N, j in 1:N
                xi = x._limbs[i]
                yj = y._limbs[j]
                p, e = MultiFloats.two_prod(xi, yj)
                @test Rational{BigInt}(p) + Rational{BigInt}(e) ==
                      Rational{BigInt}(xi) * Rational{BigInt}(yj)
            end

            z = mulN(x, y)
            terms = ExperimentalArithmetic._mul_safe_terms(x._limbs, y._limbs)
            terms_swapped = ExperimentalArithmetic._mul_safe_terms(y._limbs, x._limbs)
            full = ExperimentalArithmetic._mul_safe_full_limbs(x._limbs, y._limbs)
            exact = Rational{BigInt}(x) * Rational{BigInt}(y)
            discarded = exact_mul_tail(full, N)

            @test length(full) == 2N^2
            @test terms === terms_swapped
            @test MultiFloats.isnormalized(full)
            @test MultiFloats.isnormalized(z)
            @test Rational{BigInt}(z) + discarded == exact
            @test z === reference_mul(x, y)
            @test z === mulN(y, x)
            @test z === ExperimentalArithmetic.mul_safe(x, y)
            return z
        end

        @testset "Float64x$(N) dense/scaled" begin
            for _ in 1:35
                check_case(dense_mulN(T, N), dense_mulN(T, N))
            end
            for _ in 1:35
                check_case(
                    dense_mulN(T, N; emin=-90, emax=90),
                    dense_mulN(T, N; emin=-90, emax=90),
                )
            end
        end

        @testset "Float64x$(N) exact identities/sign" begin
            z = zero(T)
            o = one(T)
            for _ in 1:8
                x = dense_mulN(T, N; emin=-60, emax=60)
                @test mulN(x, o) === x
                @test mulN(o, x) === x
                @test iszero(mulN(x, z))
                @test iszero(mulN(z, x))

                nx = setprecision(BigFloat, 2048) do
                    T(BigFloat(-Rational{BigInt}(x)))
                end
                @test MultiFloats.isnormalized(nx)
                check_case(x, o)
                check_case(x, nx)
            end
        end

        @testset "Float64x$(N) powers and boundaries" begin
            for ex in (-100, 0, 100), ey in (-100, 0, 100)
                qx = pow2_mulN(ex)
                qy = pow2_mulN(ey)
                x = T(BigFloat(qx))
                y = T(BigFloat(qy))
                @test check_case(x, y) === T(BigFloat(qx * qy))
            end

            for e in (-60, 0, 60)
                base = pow2_mulN(e)
                for depth in (53, 53(N - 2), 53(N - 1))
                    delta = pow2_mulN(e - depth)
                    x = T(BigFloat(base + delta))
                    y = T(BigFloat(base - delta))
                    check_case(x, y)
                    check_case(x, T(BigFloat(-base + delta)))
                end
            end
        end
    end

    @testset "generic safe-mul domain" begin
        T4 = MultiFloat{Float64,4}
        T9 = MultiFloat{Float64,9}
        @test_throws ArgumentError ExperimentalArithmetic.mul_safe(T4(1.0), T4(2.0))
        @test_throws ArgumentError ExperimentalArithmetic.mul_safe(T9(1.0), T9(2.0))
    end
end
