@testset "M6 high-guard reciprocal and division reference oracles" begin
    Random.seed!(0x6d06_2026)

    ref_recip = ExperimentalArithmetic.reference_recip
    ref_div = ExperimentalArithmetic.reference_div

    function pow2_ref_nd(e::Int)
        e >= 0 ? (BigInt(1) << e) // BigInt(1) :
            BigInt(1) // (BigInt(1) << (-e))
    end

    function random_ref_nd(; words=6, emin=-300, emax=300)
        numerator = BigInt(0)
        for _ in 1:words
            numerator = (numerator << 64) + BigInt(rand(UInt64))
        end
        numerator |= BigInt(1) << (64words - 1)
        q = numerator // (BigInt(1) << 64words)
        q *= pow2_ref_nd(rand(emin:emax))
        return rand(Bool) ? q : -q
    end

    function independent_pack(::Type{T}, value::Rational{BigInt}) where {T}
        return setprecision(BigFloat, 16384) do
            T(BigFloat(value))
        end
    end

    for N in 5:8
        T = MultiFloat{Float64,N}

        @testset "Float64x$(N) scalar reciprocal/division oracle" begin
            for _ in 1:10
                x = setprecision(BigFloat, 2048) do
                    T(BigFloat(random_ref_nd()))
                end
                y = setprecision(BigFloat, 2048) do
                    T(BigFloat(random_ref_nd()))
                end
                qx = Rational{BigInt}(x)
                qy = Rational{BigInt}(y)

                rx = ref_recip(x)
                q = ref_div(x, y)

                @test rx === independent_pack(T, inv(qx))
                @test q === independent_pack(T, qx / qy)
                @test MultiFloats.isnormalized(rx)
                @test MultiFloats.isnormalized(q)
                @test rx === ref_div(one(T), x)
            end
        end

        @testset "Float64x$(N) identities and signs" begin
            o = one(T)
            z = zero(T)
            @test ref_recip(o) === o
            @test ref_div(o, o) === o

            for _ in 1:6
                x = setprecision(BigFloat, 2048) do
                    T(BigFloat(random_ref_nd(emin=-150, emax=150)))
                end
                nx = setprecision(BigFloat, 2048) do
                    T(BigFloat(-Rational{BigInt}(x)))
                end
                @test ref_div(x, o) === x
                @test ref_div(x, x) === o
                @test ref_div(z, x) === z
                # Use numerical equality for sign identities because upstream
                # unary negation may leave trailing signed-zero limbs. Oracle-vs-
                # independent-pack comparisons above remain strict bitwise gates.
                @test ref_recip(nx) == -ref_recip(x)
                @test ref_div(nx, x) == -o
                @test ref_div(x, nx) == -o
            end
        end

        @testset "Float64x$(N) exponent spread" begin
            o = one(T)
            for e in (-500, -200, 0, 200, 500)
                x = T(BigFloat(pow2_ref_nd(e)))
                rx = ref_recip(x)
                @test rx === T(BigFloat(pow2_ref_nd(-e)))
                @test ref_div(o, x) === rx
            end
        end

        @testset "Float64x$(N) Vec2 lane oracle" begin
            V = MultiFloatVec{2,Float64,N}
            for _ in 1:3
                xs = ntuple(_ -> setprecision(BigFloat, 2048) do
                    T(BigFloat(random_ref_nd(emin=-100, emax=100)))
                end, 2)
                ys = ntuple(_ -> setprecision(BigFloat, 2048) do
                    T(BigFloat(random_ref_nd(emin=-100, emax=100)))
                end, 2)
                vx = V(xs)
                vy = V(ys)
                vr = ref_recip(vx)
                vq = ref_div(vx, vy)
                for lane in 1:2
                    @test vr[lane] === ref_recip(xs[lane])
                    @test vq[lane] === ref_div(xs[lane], ys[lane])
                end
            end
        end
    end

    @testset "Float32 packing smoke" begin
        for N in (5, 8)
            T = MultiFloat{Float32,N}
            for _ in 1:4
                x = setprecision(BigFloat, 1024) do
                    T(BigFloat(random_ref_nd(words=3, emin=-40, emax=40)))
                end
                y = setprecision(BigFloat, 1024) do
                    T(BigFloat(random_ref_nd(words=3, emin=-40, emax=40)))
                end
                qx = Rational{BigInt}(x)
                qy = Rational{BigInt}(y)
                @test ref_recip(x) === independent_pack(T, inv(qx))
                @test ref_div(x, y) === independent_pack(T, qx / qy)
            end
        end
    end

    @testset "explicit reference domain" begin
        T4 = MultiFloat{Float64,4}
        T5 = MultiFloat{Float64,5}
        T9 = MultiFloat{Float64,9}
        T16 = MultiFloat{Float16,5}

        @test_throws ArgumentError ref_recip(T4(1.0))
        @test_throws ArgumentError ref_div(T9(1.0), T9(2.0))
        @test_throws ArgumentError ref_recip(T16(1.0))
        @test_throws DomainError ref_recip(T5(0.0))
        @test_throws DomainError ref_div(T5(1.0), T5(0.0))
        @test_throws DomainError ref_recip(T5(Inf))
        @test_throws DomainError ref_div(T5(Inf), T5(1.0))

        unnormalized = T5((1.0, 1.0, 0.0, 0.0, 0.0))
        @test !MultiFloats.isnormalized(unnormalized)
        @test_throws ArgumentError ref_recip(unnormalized)
        @test_throws ArgumentError ref_div(unnormalized, T5(1.0))
    end

    @testset "ambient BigFloat precision independence" begin
        T = MultiFloat{Float64,8}
        x = setprecision(BigFloat, 2048) do
            T(BigFloat(random_ref_nd(emin=-100, emax=100)))
        end
        y = setprecision(BigFloat, 2048) do
            T(BigFloat(random_ref_nd(emin=-100, emax=100)))
        end

        rlow = setprecision(BigFloat, 64) do
            ref_recip(x)
        end
        rhigh = setprecision(BigFloat, 2048) do
            ref_recip(x)
        end
        qlow = setprecision(BigFloat, 64) do
            ref_div(x, y)
        end
        qhigh = setprecision(BigFloat, 2048) do
            ref_div(x, y)
        end
        @test rlow === rhigh
        @test qlow === qhigh
    end
end
