using MultiFloatArithmetic
using MultiFloats
using Random

const E = MultiFloatArithmetic.Experimental
const T = MultiFloat{Float64,8}

function repack_inv8_bench(q::Rational{BigInt})
    return setprecision(BigFloat, 2048) do
        T(BigFloat(q))
    end
end

function moderate_inv8_bench(; emin=-80, emax=80)
    q = Rational{BigInt}(rand(T))
    e = rand(emin:emax)
    if e >= 0
        q *= BigInt(1) << e
    else
        q /= BigInt(1) << (-e)
    end
    q = rand(Bool) ? q : -q
    x = repack_inv8_bench(q)
    return iszero(x) ? T(1.25) : x
end

function residual_constant(x::T, y::T)
    # |y - 1/x| / |1/x| = |x*y - 1| exactly for rational represented values.
    residual = abs(Rational{BigInt}(x) * Rational{BigInt}(y) - one(Rational{BigInt}))
    return setprecision(BigFloat, 768) do
        BigFloat(residual) / BigFloat(2)^(-53 * 8)
    end
end

function min_time(f; samples=3)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function run_diagnostics(; n=120)
    Random.seed!(0x1a88_2026)
    xs = [moderate_inv8_bench() for _ in 1:n]

    oracle_mismatches = 0
    normalization_failures = 0
    max_seed_constant = BigFloat(0)
    max_one_constant = BigFloat(0)
    second_successes = 0
    second_domain_rejections = 0
    second_oracle_mismatches = 0

    for x in xs
        seed = E._inv8_seed(x)
        one = E.inv8_safe(x)
        oracle = E.reference_inv(x)

        oracle_mismatches += one !== oracle
        normalization_failures += !MultiFloats.isnormalized(one)
        max_seed_constant = max(max_seed_constant, residual_constant(x, seed))
        max_one_constant = max(max_one_constant, residual_constant(x, one))

        try
            two = E._inv8_safe_two_corrections(x)
            second_successes += 1
            second_oracle_mismatches += two !== oracle
        catch err
            if err isa DomainError
                second_domain_rejections += 1
            else
                rethrow()
            end
        end
    end

    println("cases=$(n)")
    println("one_correction_oracle_mismatches=$(oracle_mismatches)")
    println("one_correction_normalization_failures=$(normalization_failures)")
    println("max seed |xy-1|/u^8 = $(max_seed_constant)")
    println("max one-correction |xy-1|/u^8 = $(max_one_constant)")
    println("second_correction_successes=$(second_successes)")
    println("second_correction_domain_rejections=$(second_domain_rejections)")
    println("second_correction_oracle_mismatches=$(second_oracle_mismatches)")

    @assert oracle_mismatches == 0
    @assert normalization_failures == 0
end

function run_timing(; n=20)
    Random.seed!(0x1a80_2026)
    xs = [moderate_inv8_bench(emin=-50, emax=50) for _ in 1:n]
    out = Vector{T}(undef, n)

    candidate!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.inv8_safe(xs[i])
        end
        out
    end
    oracle!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.reference_inv(xs[i])
        end
        out
    end

    candidate!(); oracle!()
    tc = min_time(candidate!)
    to = min_time(oracle!())
    println("Float64x8 inv8_safe: $(round(tc*1e3; digits=3)) ms for $(n) cases")
    println("Float64x8 reference_inv: $(round(to*1e3; digits=3)) ms for $(n) cases")
    println("candidate/reference time ratio = $(round(tc/to; digits=3))x")
end

println("M6 Float64x8 reciprocal one-vs-two-correction diagnostic")
println("CPU: ", Sys.CPU_NAME)
run_diagnostics()
run_timing()
