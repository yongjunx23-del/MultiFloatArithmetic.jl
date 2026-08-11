using MultiFloatArithmetic
using MultiFloats
using Random
using Test

const CASES = (
    (MultiFloats.Float64x2, 2, BigFloat(34)),
    (MultiFloats.Float64x3, 3, BigFloat(184)),
    (MultiFloats.Float64x4, 4, BigFloat(812)),
)

signed_rand(::Type{T}) where {T} = rand(Bool) ? rand(T) : -rand(T)

function wide_rand(::Type{T}; emin=-400, emax=400) where {T}
    x = signed_rand(T)
    exponent = rand(emin:emax)
    return T(ldexp(big(x), exponent))
end

function check_operand_relative_bound(z, x, y, c, limbs, constant, u)
    xy = big(x) * big(y)
    reference = xy + big(c)
    scale = abs(xy) + abs(big(c))
    err = abs(big(z) - reference)
    bound = constant * u^limbs * scale
    oracle_slack = eps(BigFloat) * max(scale, one(BigFloat))
    return err <= bound + oracle_slack
end

@testset "branch-free fused FMA research kernels" begin
    Random.seed!(0x5d9a_2026)
    setprecision(BigFloat, 1024) do
        u = BigFloat(2)^(-53)
        for (T, limbs, constant) in CASES
            @testset "$(T) ordinary scalar" begin
                for _ in 1:2_000
                    x = signed_rand(T)
                    y = signed_rand(T)
                    c = signed_rand(T)
                    z = fma_fast(x, y, c)

                    @test MultiFloats.isnormalized(z)
                    @test z === fma_fast(y, x, c)
                    @test check_operand_relative_bound(z, x, y, c, limbs, constant, u)
                end
            end

            @testset "$(T) wide-exponent scalar" begin
                for _ in 1:2_000
                    x = wide_rand(T)
                    y = wide_rand(T)
                    c = wide_rand(T)
                    z = fma_fast(x, y, c)

                    @test isfinite(z)
                    @test MultiFloats.isnormalized(z)
                    @test z === fma_fast(y, x, c)
                    @test check_operand_relative_bound(z, x, y, c, limbs, constant, u)
                end
            end

            @testset "$(T) Vec4 lane equivalence" begin
                V = MultiFloatVec{4,Float64,limbs}
                for _ in 1:250
                    xs = ntuple(_ -> signed_rand(T), 4)
                    ys = ntuple(_ -> signed_rand(T), 4)
                    cs = ntuple(_ -> signed_rand(T), 4)
                    vz = fma_fast(V(xs), V(ys), V(cs))
                    for lane in 1:4
                        @test vz[lane] === fma_fast(xs[lane], ys[lane], cs[lane])
                    end
                end

                for _ in 1:250
                    xs = ntuple(_ -> wide_rand(T), 4)
                    ys = ntuple(_ -> wide_rand(T), 4)
                    cs = ntuple(_ -> wide_rand(T), 4)
                    vz = fma_fast(V(xs), V(ys), V(cs))
                    for lane in 1:4
                        scalar = fma_fast(xs[lane], ys[lane], cs[lane])
                        @test vz[lane] === scalar
                        @test check_operand_relative_bound(
                            scalar, xs[lane], ys[lane], cs[lane], limbs, constant, u)
                    end
                end
            end
        end
    end
end

@testset "SIMD width lane equivalence" begin
    Random.seed!(0x51d0_2026)
    setprecision(BigFloat, 1024) do
        u = BigFloat(2)^(-53)
        for (T, limbs, constant) in CASES
            for W in (2, 4, 8)
                V = MultiFloatVec{W,Float64,limbs}
                @testset "$(T) Vec$(W)" begin
                    for _ in 1:100
                        xs = ntuple(_ -> wide_rand(T), W)
                        ys = ntuple(_ -> wide_rand(T), W)
                        cs = ntuple(_ -> wide_rand(T), W)
                        vz = fma_fast(V(xs), V(ys), V(cs))
                        for lane in 1:W
                            scalar = fma_fast(xs[lane], ys[lane], cs[lane])
                            @test vz[lane] === scalar
                            @test check_operand_relative_bound(
                                scalar, xs[lane], ys[lane], cs[lane], limbs, constant, u)
                        end
                    end
                end
            end
        end
    end
end

@testset "destructive cancellation remains explicit" begin
    Random.seed!(0xcafe_2026)
    setprecision(BigFloat, 1024) do
        u = BigFloat(2)^(-53)
        cancellation_bits = (20, 50, 100, 150, 200)

        for (T, limbs, constant) in CASES
            x = T(BigFloat("0.812345678901234567890123456789"))
            y = T(BigFloat("0.912345678901234567890123456789"))
            c = -T(big(x) * big(y))
            z = fma_fast(x, y, c)
            @test isfinite(z)
            @test z === fma_fast(y, x, c)
            @test check_operand_relative_bound(z, x, y, c, limbs, constant, u)

            for bits in cancellation_bits
                δ = BigFloat(2)^(-bits)
                for _ in 1:100
                    x = wide_rand(T; emin=-200, emax=200)
                    y = wide_rand(T; emin=-200, emax=200)
                    xy = big(x) * big(y)
                    c = T(-xy * (one(BigFloat) - δ))
                    z = fma_fast(x, y, c)

                    @test isfinite(z)
                    @test z === fma_fast(y, x, c)
                    @test check_operand_relative_bound(z, x, y, c, limbs, constant, u)
                end
            end
        end
    end
end
