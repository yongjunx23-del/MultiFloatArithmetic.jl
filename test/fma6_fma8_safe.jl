@testset "safe Float64x6-x8 direct FMA baselines" begin
    Random.seed!(0xf668_2026)

    reference_fma = ExperimentalArithmetic.reference_fma

    function pow2_fmaN(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    function dense_fmaN(::Type{T}, N; emin=-50, emax=50) where {T}
        numerator = BigInt(0)
        for _ in 1:N
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        numerator |= BigInt(1) << (64N - 1)
        q = numerator // (BigInt(1) << 64N)
        q *= pow2_fmaN(rand(emin:emax))
        q = rand(Bool) ? q : -q
        return setprecision(BigFloat, 2048) do
            T(BigFloat(q))
        end
    end

    function exact_fma_tail(full, N)
        total = zero(Rational{BigInt})
        @inbounds for i in N+1:length(full)
            total += Rational{BigInt}(full[i])
        end
        return total
    end

    @testset "generic x5 matches frozen baseline" begin
        T5 = MultiFloat{Float64,5}
        for _ in 1:20
            x = dense_fmaN(T5, 5)
            y = dense_fmaN(T5, 5)
            c = dense_fmaN(T5, 5)
            @test ExperimentalArithmetic.fma_safe(x, y, c) ===
                  ExperimentalArithmetic.fma5_safe(x, y, c)
        end
    end

    for N in 6:8
        T = MultiFloat{Float64,N}
        fmaN = N == 6 ? ExperimentalArithmetic.fma6_safe :
               N == 7 ? ExperimentalArithmetic.fma7_safe :
                        ExperimentalArithmetic.fma8_safe

        check_case = function (x, y, c)
            @inbounds for i in 1:N, j in 1:N
                xi = x._limbs[i]
                yj = y._limbs[j]
                p, e = MultiFloats.two_prod(xi, yj)
                @test Rational{BigInt}(p) + Rational{BigInt}(e) ==
                      Rational{BigInt}(xi) * Rational{BigInt}(yj)
            end

            result = fmaN(x, y, c)
            terms = ExperimentalArithmetic._fma_safe_terms(x._limbs, y._limbs, c._limbs)
            terms_swapped = ExperimentalArithmetic._fma_safe_terms(y._limbs, x._limbs, c._limbs)
            full = ExperimentalArithmetic._fma_safe_full_limbs(x._limbs, y._limbs, c._limbs)
            exact = Rational{BigInt}(x) * Rational{BigInt}(y) + Rational{BigInt}(c)
            discarded = exact_fma_tail(full, N)

            @test length(full) == 2N^2 + N
            @test terms === terms_swapped
            @test MultiFloats.isnormalized(full)
            @test MultiFloats.isnormalized(result)
            @test Rational{BigInt}(result) + discarded == exact
            @test result === reference_fma(x, y, c)
            @test result === fmaN(y, x, c)
            @test result === ExperimentalArithmetic.fma_safe(x, y, c)
            return result
        end

        @testset "Float64x$(N) dense/scaled" begin
            ncase = N == 6 ? 24 : N == 7 ? 18 : 12
            for _ in 1:ncase
                check_case(dense_fmaN(T, N), dense_fmaN(T, N), dense_fmaN(T, N))
            end
            for _ in 1:ncase
                check_case(
                    dense_fmaN(T, N; emin=-85, emax=85),
                    dense_fmaN(T, N; emin=-85, emax=85),
                    dense_fmaN(T, N; emin=-85, emax=85),
                )
            end
        end

        @testset "Float64x$(N) identities" begin
            z = zero(T)
            o = one(T)
            for _ in 1:5
                x = dense_fmaN(T, N; emin=-60, emax=60)
                c = dense_fmaN(T, N; emin=-60, emax=60)
                @test fmaN(x, o, z) === x
                @test fmaN(x, z, c) === c
                @test fmaN(z, x, c) === c

                nx = setprecision(BigFloat, 2048) do
                    T(BigFloat(-Rational{BigInt}(x)))
                end
                @test iszero(fmaN(x, o, nx))
                check_case(x, o, nx)
            end
        end

        @testset "Float64x$(N) destructive cancellation" begin
            depths = (100, 53N - 80, 53N - 20, 53N + 30)
            ncancel = N == 6 ? 10 : N == 7 ? 8 : 6
            for depth in depths, _ in 1:ncancel
                x = dense_fmaN(T, N; emin=-40, emax=40)
                y = dense_fmaN(T, N; emin=-40, emax=40)
                qxy = Rational{BigInt}(x) * Rational{BigInt}(y)
                lead_exp = exponent(x._limbs[1]) + exponent(y._limbs[1])
                delta = pow2_fmaN(lead_exp - depth)
                c = setprecision(BigFloat, 2048) do
                    T(BigFloat(-qxy + (rand(Bool) ? delta : -delta)))
                end
                check_case(x, y, c)
            end
        end

        @testset "Float64x$(N) power boundaries" begin
            for ex in (-60, 0, 60), ey in (-60, 0, 60)
                x = T(BigFloat(pow2_fmaN(ex)))
                y = T(BigFloat(pow2_fmaN(ey)))
                c = T(BigFloat(-pow2_fmaN(ex + ey)))
                @test iszero(check_case(x, y, c))
            end
        end
    end

    @testset "generic direct-FMA domain" begin
        T4 = MultiFloat{Float64,4}
        T9 = MultiFloat{Float64,9}
        @test_throws ArgumentError ExperimentalArithmetic.fma_safe(T4(1.0), T4(2.0), T4(3.0))
        @test_throws ArgumentError ExperimentalArithmetic.fma_safe(T9(1.0), T9(2.0), T9(3.0))
    end
end
