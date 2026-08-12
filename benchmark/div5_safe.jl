using MultiFloatArithmetic
using MultiFloats
using Random
using Printf

const E = MultiFloatArithmetic.Experimental
const T5D = MultiFloat{Float64,5}

function pow2_div5_bench(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

function dense_div5_bench(; words=5, emin=-80, emax=80)
    numerator = BigInt(0)
    for _ in 1:words
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64words - 1)
    q = numerator // (BigInt(1) << 64words)
    q *= pow2_div5_bench(rand(emin:emax))
    q = rand(Bool) ? q : -q
    return setprecision(BigFloat, 1536) do
        T5D(BigFloat(q))
    end
end

function quotient_error_bits(x::T5D, y::T5D, q::T5D)
    qref = Rational{BigInt}(x) / Rational{BigInt}(y)
    err = abs(Rational{BigInt}(q) - qref)
    iszero(err) && return Inf
    return setprecision(BigFloat, 1024) do
        rel = BigFloat(err) / abs(BigFloat(qref))
        Float64(-log2(rel))
    end
end

function residual_bits_div5(x::T5D, y::T5D, q::T5D)
    r = abs(Rational{BigInt}(x) - Rational{BigInt}(q) * Rational{BigInt}(y))
    iszero(r) && return Inf
    scale = max(abs(Rational{BigInt}(x)), one(Rational{BigInt}))
    return setprecision(BigFloat, 1024) do
        Float64(-log2(BigFloat(r) / BigFloat(scale)))
    end
end

function best_ms_div5(f; samples=3)
    f()
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, 1e3 * @elapsed f())
    end
    return best
end

function main()
    Random.seed!(0x6d15_b26)
    xs = [dense_div5_bench() for _ in 1:80]
    ys = [dense_div5_bench() for _ in 1:80]

    mismatch_q0 = 0
    mismatch_q1 = 0
    min_q0_bits = Inf
    min_q1_bits = Inf
    min_res0_bits = Inf
    min_res1_bits = Inf

    for i in eachindex(xs)
        x = xs[i]
        y = ys[i]
        oracle = E.reference_div(x, y)
        q0, invy = E._div5_seed(x, y)
        q1, _ = E._div5_correct_once(x, y, q0, invy)
        mismatch_q0 += q0 !== oracle
        mismatch_q1 += q1 !== oracle
        min_q0_bits = min(min_q0_bits, quotient_error_bits(x, y, q0))
        min_q1_bits = min(min_q1_bits, quotient_error_bits(x, y, q1))
        min_res0_bits = min(min_res0_bits, residual_bits_div5(x, y, q0))
        min_res1_bits = min(min_res1_bits, residual_bits_div5(x, y, q1))
    end

    out = Vector{T5D}(undef, length(xs))
    work_q0!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E._div5_uncorrected(xs[i], ys[i])
        end
        out
    end
    work_q1!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.div5_safe(xs[i], ys[i])
        end
        out
    end

    t_q0 = best_ms_div5(work_q0!)
    t_q1 = best_ms_div5(work_q1!)

    println("M6 Float64x5 division direct-residual correction A/B")
    println("CPU: ", Sys.CPU_NAME)
    println("cases=$(length(xs))")
    println("initial x*recip5 quotient:")
    println("  oracle_mismatches=$(mismatch_q0)")
    @printf("  worst relative quotient bits ~= %.1f\n", min_q0_bits)
    @printf("  worst scaled residual bits ~= %.1f\n", min_res0_bits)
    @printf("  time = %.3f ms / %d\n", t_q0, length(xs))
    println("after one direct residual/correction:")
    println("  oracle_mismatches=$(mismatch_q1)")
    @printf("  worst relative quotient bits ~= %.1f\n", min_q1_bits)
    @printf("  worst scaled residual bits ~= %.1f\n", min_res1_bits)
    @printf("  time = %.3f ms / %d\n", t_q1, length(xs))
    @printf("corrected / uncorrected time ratio = %.3fx\n", t_q1 / t_q0)

    @assert mismatch_q1 == 0
end

main()
