@testset "specialized expansion × one-limb multiplication" begin
    Random.seed!(0x51a1_2026)

    for (T, limbs, _) in CASES
        @testset "$(T) scalar bitwise equivalence" begin
            for _ in 1:4_000
                x = wide_rand(T; emin=-300, emax=300)
                q = ldexp((rand(Bool) ? 1.0 : -1.0) * (0.5 + rand(Float64)), rand(-300:300))

                qexp = ntuple(i -> isone(i) ? q : zero(q), limbs)
                reference = MultiFloats.mfmul(x._limbs, qexp, Val(limbs))
                candidate = mul_scalar_limbs(x._limbs, q, Val(limbs))

                @test candidate === reference
            end
        end

        @testset "$(T) Vec4 bitwise equivalence" begin
            V = MultiFloatVec{4,Float64,limbs}
            V1 = MultiFloatVec{4,Float64,1}

            for _ in 1:500
                xs = ntuple(_ -> wide_rand(T; emin=-250, emax=250), 4)
                qs = ntuple(_ -> ldexp((rand(Bool) ? 1.0 : -1.0) *
                                        (0.5 + rand(Float64)), rand(-250:250)), 4)

                vx = V(xs)
                qv = V1(qs)._limbs[1]
                qexp = ntuple(i -> isone(i) ? qv : zero(qv), limbs)

                reference = MultiFloats.mfmul(vx._limbs, qexp, Val(limbs))
                candidate = mul_scalar_limbs(vx._limbs, qv, Val(limbs))

                @test candidate === reference
            end
        end
    end
end
