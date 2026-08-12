@testset "M6 reciprocal division and sqrt reference oracles" begin
    Random.seed!(0x6a61_2026)

    ref_inv = ExperimentalArithmetic.reference_inv
    ref_div = ExperimentalArithmetic.reference_div
    ref_sqrt = ExperimentalArithmetic.reference_sqrt

    function repack_test_rational(::Type{T}, q::Rational{BigInt}) where {T}
        return setprecision(BigFloat, 2048) do
            T(BigFloat(q))
        end
    end

    function moderate_value(::Type{T}; emin=-80, emax=80) where {T}
        q = Rational{BigInt}(rand(T))
        e = rand(emin:emax)
        if e >= 0
            q *= BigInt(1) << e
        else
            q /= BigInt(1) << (-e)
        end
        q = rand(Bool) ? q : -q
        return repack_test_rational(T, q)
    end

    function positive_value(::Type{T}; emin=-80, emax=80) where {T}
        q = abs(Rational{BigInt}(moderate_value(T; emin=emin, emax=emax)))
        return repack_test_rational(T, q)
    end

    function hp_inv(::Type{T}, x; p=16384) where {T}
        q = Rational{BigInt}(x)
        return setprecision(BigFloat, p) do
            T(BigFloat(denominator(q)) / BigFloat(numerator(q)))
        end
    end

    function hp_div(::Type{T}, x, y; p=16384) where {T}
        q = Rational{BigInt}(x) / Rational{BigInt}(y)
        return setprecision(BigFloat, p) do
            T(BigFloat(numerator(q)) / BigFloat(denominator(q)))
        end
    end

    function hp_sqrt(::Type{T}, x; p=16384) where {T}
        q = Rational{BigInt}(x)
        return setprecision(BigFloat, p) do
            T(sqrt(BigFloat(numerator(q)) / BigFloat(denominator(q))))
        end
    end

    for N in 5:8
        T = MultiFloat{Float64,N}
        @testset "Float64x$(N) scalar oracle" begin
            for _ in 1:10
                x = moderate_value(T)
                y = moderate_value(T)
                iszero(x) && (x = T(1.25))
                iszero(y) && (y = T(0.75))
                positive = positive_value(T)
                iszero(positive) && (positive = T(2.0))

                invx = ref_inv(x)
                quot = ref_div(x, y)
                root = ref_sqrt(positive)

                @test invx === hp_inv(T, x)
                @test quot === hp_div(T, x, y)
                @test root === hp_sqrt(T, positive)
                @test MultiFloats.isnormalized(invx)
                @test MultiFloats.isnormalized(quot)
                @test MultiFloats.isnormalized(root)
            end
        end

        @testset "Float64x$(N) exact identities" begin
            z = zero(T)
            o = one(T)
            two = T(2.0)
            four = T(4.0)
            x = moderate_value(T; emin=-40, emax=40)
            iszero(x) && (x = T(1.25))

            @test ref_inv(o) === o
            @test ref_div(x, o) === x
            @test ref_div(x, x) === o
            @test ref_sqrt(z) === z
            @test ref_sqrt(o) === o
            @test ref_sqrt(four) === two
        end

        @testset "Float64x$(N) Vec2 lane oracle" begin
            V = MultiFloatVec{2,Float64,N}
            xs = ntuple(_ -> begin
                x = moderate_value(T; emin=-40, emax=40)
                iszero(x) ? T(1.25) : x
            end, 2)
            ys = ntuple(_ -> begin
                y = moderate_value(T; emin=-40, emax=40)
                iszero(y) ? T(0.75) : y
            end, 2)
            ps = ntuple(_ -> begin
                x = positive_value(T; emin=-40, emax=40)
                iszero(x) ? T(2.0) : x
            end, 2)

            vi = ref_inv(V(xs))
            vd = ref_div(V(xs), V(ys))
            vs = ref_sqrt(V(ps))
            for lane in 1:2
                @test vi[lane] === ref_inv(xs[lane])
                @test vd[lane] === ref_div(xs[lane], ys[lane])
                @test vs[lane] === ref_sqrt(ps[lane])
            end
        end
    end

    @testset "Float32 x5/x8 smoke" begin
        for N in (5, 8)
            T = MultiFloat{Float32,N}
            x = T(Float32(1.234567))
            y = T(Float32(0.7654321))
            p = T(Float32(2.345678))
            @test ref_inv(x) === hp_inv(T, x)
            @test ref_div(x, y) === hp_div(T, x, y)
            @test ref_sqrt(p) === hp_sqrt(T, p)
        end
    end

    @testset "domain and ambient precision" begin
        T5 = MultiFloat{Float64,5}
        T4 = MultiFloat{Float64,4}
        T9 = MultiFloat{Float64,9}

        @test_throws DomainError ref_inv(zero(T5))
        @test_throws DomainError ref_div(one(T5), zero(T5))
        @test_throws DomainError ref_sqrt(T5(-1.0))
        @test_throws DomainError ref_inv(T5(Inf))
        @test_throws ArgumentError ref_inv(T4(1.0))
        @test_throws ArgumentError ref_sqrt(T9(1.0))

        x = T5(BigFloat("1.234567890123456789"))
        low = setprecision(BigFloat, 64) do
            (ref_inv(x), ref_div(x, T5(0.75)), ref_sqrt(x))
        end
        high = setprecision(BigFloat, 4096) do
            (ref_inv(x), ref_div(x, T5(0.75)), ref_sqrt(x))
        end
        @test low === high
    end
end
