function _paper_p_nonoverlap(limbs::NTuple{3,Float64})
    @inbounds for i in 1:2
        hi = limbs[i]
        lo = limbs[i + 1]
        if iszero(hi)
            iszero(lo) || return false
        elseif !(abs(lo) < eps(abs(hi)))
            return false
        end
    end
    return true
end

function _relative_product_error(z, x, y)
    reference = big(x) * big(y)
    iszero(reference) && return zero(BigFloat)
    return abs(big(z) - reference) / abs(reference)
end

@testset "Fabiano-Muller-Picot 2019 fast triple-word products" begin
    Random.seed!(0x3a19_2019)
    setprecision(BigFloat, 1024) do
        u = BigFloat(2)^(-53)
        bound33 = 44 * u^3 + 176 * u^4
        bound23 = 18 * u^3 + 75 * u^4
        oracle_slack = BigFloat(2)^(-900)

        @testset "3Prod_fast(3,3) paper bound" begin
            for _ in 1:20_000
                x = rand(Float64x3) + Float64x3(0.5)
                y = rand(Float64x3) + Float64x3(0.5)
                z = tw_prod33_fast(x, y)

                @test _relative_product_error(z, x, y) <= bound33 + oracle_slack
                @test _paper_p_nonoverlap(z._limbs)
                @test z === tw_prod33_fast(y, x)
            end

            for _ in 1:5_000
                x = wide_rand(Float64x3; emin=-200, emax=200)
                y = wide_rand(Float64x3; emin=-200, emax=200)
                iszero(x) && continue
                iszero(y) && continue
                z = tw_prod33_fast(x, y)

                @test isfinite(z)
                @test _relative_product_error(z, x, y) <= bound33 + oracle_slack
                @test _paper_p_nonoverlap(z._limbs)
                @test z === tw_prod33_fast(y, x)
            end
        end

        @testset "3Prod_fast(2,3) paper bound" begin
            for _ in 1:20_000
                x = rand(Float64x2) + Float64x2(0.5)
                y = rand(Float64x3) + Float64x3(0.5)
                z = tw_prod23_fast(x, y)

                @test _relative_product_error(z, x, y) <= bound23 + oracle_slack
                @test _paper_p_nonoverlap(z._limbs)
            end

            for _ in 1:5_000
                x = wide_rand(Float64x2; emin=-200, emax=200)
                y = wide_rand(Float64x3; emin=-200, emax=200)
                iszero(x) && continue
                iszero(y) && continue
                z = tw_prod23_fast(x, y)

                @test isfinite(z)
                @test _relative_product_error(z, x, y) <= bound23 + oracle_slack
                @test _paper_p_nonoverlap(z._limbs)
            end
        end

        @testset "MultiFloats strong-normalization compatibility probe" begin
            # This is a hard compatibility gate for using a literature TW result
            # as an ordinary MultiFloat.  The 2019 paper only proves the weaker
            # P-nonoverlap invariant; if this test fails we must keep the TW
            # kernels behind a separate representation or explicitly renormalize.
            for _ in 1:10_000
                x = rand(Float64x3) + Float64x3(0.5)
                y = rand(Float64x3) + Float64x3(0.5)
                @test MultiFloats.isnormalized(tw_prod33_fast(x, y))
            end
        end
    end
end
