using MultiFloatArithmetic
using MultiFloats
using Random
using Printf

const E = MultiFloatArithmetic.Experimental
const T5R = MultiFloat{Float64,5}

function pow2_recip5_bench(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

function dense_recip5_bench(; words=5, emin=-80, emax=80)
    numerator = BigInt(0)
    for _ in 1:words
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64words - 1)
    q = numerator // (BigInt(1) << 64words)
    q *= pow2_recip5_bench(rand(emin:emax))
    q = rand(Bool) ? q : -q
    return setprecision(BigFloat, 1536) do
        T5R(BigFloat(q))
    end
end

function residual_bits(x::T5R, z::T5R)
    r = abs(one(Rational{BigInt}) - Rational{BigInt}(x) * Rational{BigInt}(z))
    iszero(r) && return Inf
    return setprecision(BigFloat, 1024) do
        Float64(-log2(BigFloat(r)))
    end
end

function best_ms(f; samples=3)
    f()
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, 1e3 * @elapsed f())
    end
    return best
end

Random.seed!(0x6ec5_b26)
xs = [dense_recip5_bench() for _ in 1:120]

mismatch_x4 = 0
mismatch_f64 = 0
min_seed_bits = Inf
min_x4one_bits = Inf
min_f64three_bits = Inf
for x in xs
    oracle = E.reference_inv(x)
    seed = E._recip5_seed_x4(x)
    x4one = E._recip5_x4_one_correction(x)
    f64three = E._recip5_float64_three_corrections(x)
    mismatch_x4 += x4one !== oracle
    mismatch_f64 += f64three !== oracle
    min_seed_bits = min(min_seed_bits, residual_bits(x, seed))
    min_x4one_bits = min(min_x4one_bits, residual_bits(x, x4one))
    min_f64three_bits = min(min_f64three_bits, residual_bits(x, f64three))
end

out = Vector{T5R}(undef, length(xs))
work_x4!() = begin
    @inbounds for i in eachindex(xs)
        out[i] = E._recip5_x4_one_correction(xs[i])
    end
    out
end
work_f64!() = begin
    @inbounds for i in eachindex(xs)
        out[i] = E._recip5_float64_three_corrections(xs[i])
    end
    out
end
work_public!() = begin
    @inbounds for i in eachindex(xs)
        out[i] = E.recip5_safe(xs[i])
    end
    out
end

t_x4 = best_ms(work_x4!)
t_f64 = best_ms(work_f64!)
t_public = best_ms(work_public!)

println("M6 Float64x5 reciprocal seed/correction A/B")
println("CPU: ", Sys.CPU_NAME)
println("cases=$(length(xs))")
println("x4 seed + 1 x5 direct correction:")
println("  oracle_mismatches=$(mismatch_x4)")
@printf("  worst exact residual bits ~= %.1f\n", min_x4one_bits)
@printf("  time = %.3f ms / %d\n", t_x4, length(xs))
println("Float64 seed + 3 x5 direct corrections:")
println("  oracle_mismatches=$(mismatch_f64)")
@printf("  worst exact residual bits ~= %.1f\n", min_f64three_bits)
@printf("  time = %.3f ms / %d\n", t_f64, length(xs))
@printf("x4 seed before x5 correction worst residual bits ~= %.1f\n", min_seed_bits)
@printf("three-correction / one-correction time ratio = %.3fx\n", t_f64 / t_x4)
@printf("public recip5_safe time = %.3f ms / %d\n", t_public, length(xs))

@assert mismatch_x4 == 0
@assert mismatch_f64 == 0
@assert min_x4one_bits >= 250
@assert min_f64three_bits >= 250
