using MultiFloatArithmetic
using MultiFloats
using Random
using Test

include(joinpath(@__DIR__, "..", "audit", "paper2607_v1.jl"))
using .Paper2607V1Audit

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
    truth = rx * ry + rc
    scale = abs(rx * ry) + abs(rc)
    err = abs(paper_exact(z) - truth)
    iszero(scale) && return iszero(err)
    return setprecision(BigFloat, 1024) do
        u = BigFloat(2)^(-53)
        BigFloat(err) <= BigFloat(C) * u^K * BigFloat(scale)
    end
end

@testset "arXiv:2607.11391v1 DW Algorithm 1" begin
    Random.seed!(0x2607_0002)
    T = Float64x2
    for _ in 1:2_000
        x = paper_wide_rand(T)
        y = paper_wide_rand(T)
        c = paper_wide_rand(T)
        paper = paper_dw_fma(x, y, c)
        @test paper === fma_fast(x, y, c)
        @test MultiFloats.isnormalized(paper)
        @test paper === paper_dw_fma(y, x, c)
        @test paper_bound_ok(paper, x, y, c, 2, 34)
    end
end

@testset "arXiv:2607.11391v1 QW Algorithm 3 ordinary/wide" begin
    Random.seed!(0x2607_0004)
    T = Float64x4
    for generator in (paper_signed_rand, paper_wide_rand)
        for _ in 1:1_000
            x = generator(T)
            y = generator(T)
            c = generator(T)
            paper = paper_qw_fma(x, y, c)
            safe = fma_fast(x, y, c)

            @test all(paper_qw_fast_preconditions(x, y, c))
            @test paper === paper_qw_fma(y, x, c)
            @test paper_exact(paper) == paper_exact(safe)
            @test paper_bound_ok(paper, x, y, c, 4, 812)
        end
    end
end

@testset "arXiv:2607.11391v1 QW cancellation non-overlap caveat" begin
    Random.seed!(0x2607_cafe)
    T = Float64x4
    nonnormalized = 0
    changed_by_canonicalization = 0
    cases = 0
    setprecision(BigFloat, 1024) do
        for bits in (20, 50, 80, 110, 150)
            delta = BigFloat(2)^(-bits)
            for _ in 1:200
                x = paper_wide_rand(T; emin=-180, emax=180)
                y = paper_wide_rand(T; emin=-180, emax=180)
                xy = BigFloat(x) * BigFloat(y)
                c = T(-xy * (one(BigFloat) - delta))
                paper = paper_qw_fma(x, y, c)
                safe = fma_fast(x, y, c)

                @test all(paper_qw_fast_preconditions(x, y, c))
                @test paper === paper_qw_fma(y, x, c)
                @test paper_exact(paper) == paper_exact(safe)
                @test paper_bound_ok(paper, x, y, c, 4, 812)

                nonnormalized += !MultiFloats.isnormalized(paper)
                changed_by_canonicalization += (paper !== safe)
                cases += 1
            end
        end
    end

    # This is the paper's Sec. 2.9.5 caveat made concrete: the cheap QW value
    # remains accurate, but its four-component representation need not satisfy
    # a consumer's strict non-overlap/canonical contract after cancellation.
    @test nonnormalized > 0
    @test changed_by_canonicalization > 0
    @test cases == 1_000
end
