@testset "x5-x8 exact-rational reference arithmetic" begin
    Random.seed!(0x58a2_2026)

    function random_dyadic(; words=8, emin=-400, emax=400)
        numerator = BigInt(0)
        for _ in 1:words
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        bits = 64 * words
        numerator = numerator | (BigInt(1) << (bits - 1))
        value = numerator // (BigInt(1) << bits)
        exponent = rand(emin:emax)
        if exponent >= 0
            value *= BigInt(1) << exponent
        else
            value /= BigInt(1) << (-exponent)
        end
        return rand(Bool) ? value : -value
    end

    function high_precision_pack(::Type{T}, value::Rational{BigInt}) where {T}
        return setprecision(BigFloat, 8192) do
            T(BigFloat(value))
        end
    end

    ref_add = ExperimentalArithmetic.reference_add
    ref_sub = ExperimentalArithmetic.reference_sub
    ref_mul = ExperimentalArithmetic.reference_mul
    ref_fma = ExperimentalArithmetic.reference_fma

    for N in 5:8
        T = MultiFloat{Float64,N}
        @testset "Float64x$(N) scalar oracle" begin
            for _ in 1:12
                x = T(random_dyadic())
                y = T(random_dyadic())
                c = T(random_dyadic())

                qx = Rational{BigInt}(x)
                qy = Rational{BigInt}(y)
                qc = Rational{BigInt}(c)

                add = ref_add(x, y)
                sub = ref_sub(x, y)
                mul = ref_mul(x, y)
                fused = ref_fma(x, y, c)

                @test add === high_precision_pack(T, qx + qy)
                @test sub === high_precision_pack(T, qx - qy)
                @test mul === high_precision_pack(T, qx * qy)
                @test fused === high_precision_pack(T, qx * qy + qc)

                @test MultiFloats.isnormalized(add)
                @test MultiFloats.isnormalized(sub)
                @test MultiFloats.isnormalized(mul)
                @test MultiFloats.isnormalized(fused)
                @test mul === ref_mul(y, x)
                @test fused === ref_fma(y, x, c)
            end
        end

        @testset "Float64x$(N) exact identities and cancellation" begin
            x = T(random_dyadic(emin=-100, emax=100))
            z = zero(T)
            o = one(T)

            @test ref_add(x, -x) === z
            @test ref_sub(x, x) === z
            @test ref_mul(x, o) === x
            @test ref_fma(x, o, -x) === z

            y = T(random_dyadic(emin=-100, emax=100))
            exact_product = Rational{BigInt}(x) * Rational{BigInt}(y)
            c = T(-exact_product)
            fused = ref_fma(x, y, c)
            expected = high_precision_pack(
                T,
                exact_product + Rational{BigInt}(c),
            )
            @test fused === expected
            @test MultiFloats.isnormalized(fused)
        end

        @testset "Float64x$(N) Vec2 lane oracle" begin
            W = 2
            V = MultiFloatVec{W,Float64,N}
            for _ in 1:4
                xs = ntuple(_ -> T(random_dyadic(emin=-200, emax=200)), W)
                ys = ntuple(_ -> T(random_dyadic(emin=-200, emax=200)), W)
                cs = ntuple(_ -> T(random_dyadic(emin=-200, emax=200)), W)
                vx = V(xs)
                vy = V(ys)
                vc = V(cs)

                va = ref_add(vx, vy)
                vm = ref_mul(vx, vy)
                vf = ref_fma(vx, vy, vc)
                for lane in 1:W
                    @test va[lane] === ref_add(xs[lane], ys[lane])
                    @test vm[lane] === ref_mul(xs[lane], ys[lane])
                    @test vf[lane] === ref_fma(xs[lane], ys[lane], cs[lane])
                end
            end
        end
    end

    @testset "reference domain is explicit" begin
        T4 = MultiFloat{Float64,4}
        T9 = MultiFloat{Float64,9}
        @test_throws ArgumentError ref_add(T4(1.0), T4(2.0))
        @test_throws ArgumentError ref_mul(T9(1.0), T9(2.0))

        T5 = MultiFloat{Float64,5}
        @test_throws DomainError ref_add(T5(Inf), T5(1.0))

        unnormalized = T5((1.0, 1.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError ref_add(unnormalized, T5(1.0))
    end

    @testset "reference arithmetic ignores ambient BigFloat precision" begin
        T = MultiFloat{Float64,8}
        x = T(random_dyadic(emin=-100, emax=100))
        y = T(random_dyadic(emin=-100, emax=100))
        c = T(random_dyadic(emin=-100, emax=100))

        low = setprecision(BigFloat, 64) do
            ref_fma(x, y, c)
        end
        high = setprecision(BigFloat, 2048) do
            ref_fma(x, y, c)
        end
        @test low === high
    end
end
