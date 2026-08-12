using MultiFloatArithmetic
using MultiFloats
using Random
using Test

const ExperimentalArithmetic = MultiFloatArithmetic.Experimental

const CASES = (
    (MultiFloats.Float64x2, 2, BigFloat(34)),
    (MultiFloats.Float64x3, 3, BigFloat(184)),
    (MultiFloats.Float64x4, 4, BigFloat(812)),
)

const FLOAT32_CASES = (
    (MultiFloats.Float32x2, 2, BigFloat(34)),
    (MultiFloats.Float32x3, 3, BigFloat(184)),
    (MultiFloats.Float32x4, 4, BigFloat(812)),
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

@testset "public API boundary" begin
    public_names = names(MultiFloatArithmetic; all=false, imported=false)
    experimental_names = names(ExperimentalArithmetic; all=false, imported=false)

    @test :fma_fast in public_names
    @test :fma_fast_limbs in public_names
    @test :Experimental in public_names
    @test :fma5_safe ∉ public_names
    @test :add_safe ∉ public_names
    @test :add5_safe ∉ public_names
    @test :add6_safe ∉ public_names
    @test :add7_safe ∉ public_names
    @test :add8_safe ∉ public_names
    @test :mul_safe ∉ public_names
    @test :mul5_safe ∉ public_names
    @test :mul6_safe ∉ public_names
    @test :mul7_safe ∉ public_names
    @test :mul8_safe ∉ public_names
    @test :div_digits ∉ public_names
    @test :mul_scalar ∉ public_names
    @test :reference_add ∉ public_names
    @test :reference_recip ∉ public_names
    @test :reference_div ∉ public_names
    @test :fma5_safe in experimental_names
    @test :add_safe in experimental_names
    @test :add5_safe in experimental_names
    @test :add6_safe in experimental_names
    @test :add7_safe in experimental_names
    @test :add8_safe in experimental_names
    @test :mul_safe in experimental_names
    @test :mul5_safe in experimental_names
    @test :mul6_safe in experimental_names
    @test :mul7_safe in experimental_names
    @test :mul8_safe in experimental_names
    @test :div_digits in experimental_names
    @test :mul_scalar in experimental_names
    @test :reference_add in experimental_names
    @test :reference_sub in experimental_names
    @test :reference_mul in experimental_names
    @test :reference_fma in experimental_names
    @test :reference_recip in experimental_names
    @test :reference_div in experimental_names

    T5 = MultiFloat{Float64,5}
    @test_throws ArgumentError fma_fast(T5(1.0), T5(1.0), T5(0.0))
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

            @testset "$(T) exact identities" begin
                z = zero(T)
                o = one(T)
                @test fma_fast(z, o, o) == o
                @test fma_fast(o, o, z) == o
                @test iszero(fma_fast(o, o, -o))
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

@testset "Float32 implementation smoke coverage" begin
    Random.seed!(0xf032_2026)
    setprecision(BigFloat, 512) do
        u = BigFloat(2)^(-24)
        for (T, limbs, constant) in FLOAT32_CASES
            for _ in 1:500
                x = wide_rand(T; emin=-40, emax=40)
                y = wide_rand(T; emin=-40, emax=40)
                c = wide_rand(T; emin=-40, emax=40)
                z = fma_fast(x, y, c)

                @test isfinite(z)
                @test MultiFloats.isnormalized(z)
                @test z === fma_fast(y, x, c)
                @test check_operand_relative_bound(z, x, y, c, limbs, constant, u)
            end

            V = MultiFloatVec{4,Float32,limbs}
            for _ in 1:100
                xs = ntuple(_ -> wide_rand(T; emin=-30, emax=30), 4)
                ys = ntuple(_ -> wide_rand(T; emin=-30, emax=30), 4)
                cs = ntuple(_ -> wide_rand(T; emin=-30, emax=30), 4)
                vz = fma_fast(V(xs), V(ys), V(cs))
                for lane in 1:4
                    @test vz[lane] === fma_fast(xs[lane], ys[lane], cs[lane])
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

include("mul_scalar.jl")
include("division_digits.jl")
include("reference_arithmetic.jl")
include("reference_nondyadic.jl")
include("add5_safe.jl")
include("add6_add8_safe.jl")
include("mul5_safe.jl")
include("mul6_mul8_safe.jl")
include("fma5_safe.jl")
include("paper2607_v4.jl")
