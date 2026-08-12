using MultiFloatArithmetic
using MultiFloats
using Random
using Printf

const E = MultiFloatArithmetic.Experimental
const T8I = MultiFloat{Float64,8}

function pow2_inv8_bench(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

function dense_inv8_bench(; words=8, emin=-80, emax=80)
    numerator = BigInt(0)
    for _ in 1:words
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64words - 1)
    q = numerator // (BigInt(1) << 64words)
    q *= pow2_inv8_bench(rand(emin:emax))
    q = rand(Bool) ? q : -q
    return setprecision(BigFloat, 2048) do
        T8I(BigFloat(q))
    end
end

function residual_bits_inv8(x::T8I, y::T8I)
    r = abs(one(Rational{BigInt}) - Rational{BigInt}(x) * Rational{BigInt}(y))
    iszero(r) && return Inf
    return setprecision(BigFloat, 1536) do
        Float64(-log2(BigFloat(r)))
    end
end

function best_ms_inv8(f; samples=3)
    f()
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, 1e3 * @elapsed f())
    end
    return best
end

function main()
    Random.seed!(0x6ec8_b26)
    xs = [dense_inv8_bench() for _ in 1:16]

    mismatch_one = 0
    mismatch_two = 0
    min_seed_bits = Inf
    min_one_bits = Inf
    min_two_bits = Inf
    for x in xs
        oracle = E.reference_inv(x)
        seed = E._inv8_seed_x4(x)
        one_step = E._inv8_direct_one_correction(x)
        two_step = E._inv8_direct_two_corrections(x)
        mismatch_one += one_step !== oracle
        mismatch_two += two_step !== oracle
        min_seed_bits = min(min_seed_bits, residual_bits_inv8(x, seed))
        min_one_bits = min(min_one_bits, residual_bits_inv8(x, one_step))
        min_two_bits = min(min_two_bits, residual_bits_inv8(x, two_step))
    end

    out = Vector{T8I}(undef, length(xs))
    work_one!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E._inv8_direct_one_correction(xs[i])
        end
        out
    end
    work_two!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E._inv8_direct_two_corrections(xs[i])
        end
        out
    end
    work_public!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.inv8_safe(xs[i])
        end
        out
    end

    t_one = best_ms_inv8(work_one!)
    t_two = best_ms_inv8(work_two!)
    t_public = best_ms_inv8(work_public!)

    println("M6 Float64x8 direct-FMA reciprocal correction-count A/B")
    println("CPU: ", Sys.CPU_NAME)
    println("cases=$(length(xs))")
    println("x4 seed + 1 direct x8 correction:")
    println("  oracle_mismatches=$(mismatch_one)")
    @printf("  worst exact residual bits ~= %.1f\n", min_one_bits)
    @printf("  time = %.3f ms / %d\n", t_one, length(xs))
    println("x4 seed + 2 direct x8 corrections:")
    println("  oracle_mismatches=$(mismatch_two)")
    @printf("  worst exact residual bits ~= %.1f\n", min_two_bits)
    @printf("  time = %.3f ms / %d\n", t_two, length(xs))
    @printf("x4 seed before x8 correction worst residual bits ~= %.1f\n", min_seed_bits)
    @printf("two-correction / one-correction time ratio = %.3fx\n", t_two / t_one)
    @printf("public inv8_safe time = %.3f ms / %d\n", t_public, length(xs))

    @assert mismatch_two == 0
    @assert min_two_bits >= 410
end

main()
