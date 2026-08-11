@testset "quotient-digit division research candidate" begin
    Random.seed!(0xd1a1_2026)
    setprecision(BigFloat, 1024) do
        for (T, limbs, _) in CASES
            target_bits = 53 * limbs - 8
            tolerance = BigFloat(2)^(-target_bits)

            @testset "$(T) scalar accuracy" begin
                for _ in 1:2_000
                    ex = rand(-300:300)
                    ey = rand(-300:300)
                    mx = BigFloat(0.5 + rand(Float64))
                    my = BigFloat(0.5 + rand(Float64))
                    x = T(ldexp(rand(Bool) ? mx : -mx, ex))
                    y = T(ldexp(rand(Bool) ? my : -my, ey))

                    z = div_digits(x, y)
                    reference = big(x) / big(y)
                    relative_error = abs(big(z) - reference) / abs(reference)

                    @test isfinite(z)
                    @test MultiFloats.isnormalized(z)
                    @test relative_error <= tolerance
                end
            end

            @testset "$(T) Vec4 lane equivalence" begin
                V = MultiFloatVec{4,Float64,limbs}
                for _ in 1:250
                    xs = ntuple(4) do _
                        e = rand(-200:200)
                        m = BigFloat(0.5 + rand(Float64))
                        T(ldexp(rand(Bool) ? m : -m, e))
                    end
                    ys = ntuple(4) do _
                        e = rand(-200:200)
                        m = BigFloat(0.5 + rand(Float64))
                        T(ldexp(rand(Bool) ? m : -m, e))
                    end

                    vz = div_digits(V(xs), V(ys))
                    for lane in 1:4
                        scalar = div_digits(xs[lane], ys[lane])
                        @test vz[lane] === scalar
                        reference = big(xs[lane]) / big(ys[lane])
                        relative_error = abs(big(scalar) - reference) / abs(reference)
                        @test relative_error <= tolerance
                    end
                end
            end
        end
    end
end
