using MultiFloatArithmetic
using MultiFloats
using Random
using Test

include(joinpath(@__DIR__, "..", "audit", "paper2607_v4.jl"))
using .Paper2607Audit

paper_exact(x) = Rational{BigInt}(x)

function paper_signed_rand(::Type{T}) where {T}
    x = rand(T)
    return MultiFloats.renormalize(rand(Bool) ? x : -x)
end

function paper_wide_rand(::Type{T}; emin=-350, emax=350) where {T}
    return MultiFloats.renormalize(ldexp(paper_signed_rand(T), rand(emin:emax)))
end

function paper_bound_ok(z, x, y, c, K, C)
    rx = paper_exact(x)
    ry = paper_exact(y)
    rc = paper_exact(c)
    rz = paper_exact(z)
    truth = rx * ry + rc
    scale = abs(rx * ry) + abs(rc)
    err = abs(rz - truth)
    iszero(scale) && return iszero(err)
    return setprecision(BigFloat, 1024) do
        u = BigFloat(2)^(-53)
        BigFloat(err) <= BigFloat(C) * u^K * BigFloat(scale)
    end
end

@testset "arXiv:2607.11391v4 DW reproduction" begin
    Random.seed!(0x2607_0002)
    T = Float64x2
    for _ in 1:2_000
        x = paper_wide_rand(T)
        y = paper_wide_rand(T)
        c = paper_wide_rand(T)
        paper = paper_dw_v4(x, y, c)
        current = fma_fast(x, y, c)
        @test paper === current
        @test MultiFloats.isnormalized(paper)
        @test paper_bound_ok(paper, x, y, c, 2, 35)
    end
end

@testset "arXiv:2607.11391v4 QW ordinary/wide reproduction" begin
    Random.seed!(0x2607_0004)
    T = Float64x4
    for generator in (paper_signed_rand, paper_wide_rand)
        for _ in 1:1_000
            x = generator(T)
            y = generator(T)
            c = generator(T)
            v1 = paper_qw_v1_2pass(x, y, c)
            v4 = paper_qw_v4_5pass(x, y, c)
            safe = fma_fast(x, y, c)

            @test all(paper_qw_v4_fast_preconditions(x, y, c))
            @test MultiFloats.isnormalized(v4)
            @test v4 === paper_qw_v4_5pass(y, x, c)
            @test paper_exact(v4) == paper_exact(safe)
            @test paper_exact(v1) == paper_exact(v4)
            @test paper_bound_ok(v4, x, y, c, 4, 822)
        end
    end
end

@testset "arXiv:2607.11391v4 QW destructive cancellation" begin
    Random.seed!(0x2607_cafe)
    T = Float64x4
    setprecision(BigFloat, 1024) do
        for bits in (20, 50, 80, 110)
            delta = BigFloat(2)^(-bits)
            for _ in 1:250
                x = paper_wide_rand(T; emin=-180, emax=180)
                y = paper_wide_rand(T; emin=-180, emax=180)
                xy = BigFloat(x) * BigFloat(y)
                c = T(-xy * (one(BigFloat) - delta))
                v1 = paper_qw_v1_2pass(x, y, c)
                v4 = paper_qw_v4_5pass(x, y, c)
                safe = fma_fast(x, y, c)

                @test all(paper_qw_v4_fast_preconditions(x, y, c))
                @test MultiFloats.isnormalized(v4)
                @test paper_exact(v4) == paper_exact(safe)
                @test paper_exact(v1) == paper_exact(v4)
                @test paper_bound_ok(v4, x, y, c, 4, 822)
            end
        end
    end
end
