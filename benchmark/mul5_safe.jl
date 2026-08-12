using MultiFloatArithmetic
using MultiFloats
using Random

const E = MultiFloatArithmetic.Experimental
const T = MultiFloat{Float64,5}

function pow2_mul5_bench(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

function dense_mul5_bench(; words=5, emin=-100, emax=100)
    numerator = BigInt(0)
    for _ in 1:words
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64words - 1)
    q = numerator // (BigInt(1) << 64words)
    q *= pow2_mul5_bench(rand(emin:emax))
    q = rand(Bool) ? q : -q
    return setprecision(BigFloat, 1024) do
        T(BigFloat(q))
    end
end

mutable struct Mul5Stats
    n::Int
    mismatches::Int
    normalization_failures::Int
    commutativity_failures::Int
    nonzero_tail::Int
    max_relative_constant::BigFloat
end
Mul5Stats() = Mul5Stats(0, 0, 0, 0, 0, BigFloat(0))

function evaluate!(stats::Mul5Stats, x::T, y::T)
    z = E.mul5_safe(x, y)
    oracle = E.reference_mul(x, y)
    qx = Rational{BigInt}(x)
    qy = Rational{BigInt}(y)
    qz = Rational{BigInt}(z)
    exact = qx * qy
    err = abs(qz - exact)

    stats.n += 1
    stats.mismatches += z !== oracle
    stats.normalization_failures += !MultiFloats.isnormalized(z)
    stats.commutativity_failures += z !== E.mul5_safe(y, x)
    stats.nonzero_tail += !iszero(err)

    if !iszero(err) && !iszero(exact)
        setprecision(BigFloat, 768) do
            constant = BigFloat(err) / (BigFloat(2)^(-53 * 5) * abs(BigFloat(exact)))
            stats.max_relative_constant = max(stats.max_relative_constant, constant)
        end
    end
    return nothing
end

function report(label, stats::Mul5Stats)
    println(label)
    println("  cases=$(stats.n)")
    println("  oracle_mismatches=$(stats.mismatches)")
    println("  normalization_failures=$(stats.normalization_failures)")
    println("  commutativity_failures=$(stats.commutativity_failures)")
    println("  nonzero_discarded_tail=$(stats.nonzero_tail)")
    println("  max |err|/(u^5*|x*y|) = $(stats.max_relative_constant)")
end

function assert_gate(stats::Mul5Stats)
    @assert stats.mismatches == 0
    @assert stats.normalization_failures == 0
    @assert stats.commutativity_failures == 0
    # Initial empirical gate only; M4 formal analysis must replace or justify it.
    @assert stats.max_relative_constant <= 1
end

function min_time(f; samples=5)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function run_diagnostics()
    Random.seed!(0x4d55_2026)

    ordinary = Mul5Stats()
    for _ in 1:300
        evaluate!(ordinary, dense_mul5_bench(emin=-40, emax=40), dense_mul5_bench(emin=-40, emax=40))
    end
    report("dense ordinary", ordinary)
    assert_gate(ordinary)

    scaled = Mul5Stats()
    for _ in 1:300
        evaluate!(scaled, dense_mul5_bench(emin=-120, emax=120), dense_mul5_bench(emin=-120, emax=120))
    end
    report("scaled", scaled)
    assert_gate(scaled)

    boundaries = Mul5Stats()
    for ex in (-120, -40, 0, 40, 120), ey in (-120, -40, 0, 40, 120)
        x = T(BigFloat(pow2_mul5_bench(ex)))
        y = T(BigFloat(pow2_mul5_bench(ey)))
        evaluate!(boundaries, x, y)
    end
    report("power-of-two boundaries", boundaries)
    assert_gate(boundaries)
end

function run_timing(; n=500)
    Random.seed!(0x4d50_2026)
    xs = [dense_mul5_bench(emin=-40, emax=40) for _ in 1:n]
    ys = [dense_mul5_bench(emin=-40, emax=40) for _ in 1:n]
    out = Vector{T}(undef, n)

    work!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.mul5_safe(xs[i], ys[i])
        end
        out
    end

    work!()
    t = min_time(work!)
    println("Float64x5 safe multiplication: $(round(t*1e3; digits=3)) ms for $(n) cases")
end

println("M4 Float64x5 safe-multiplication diagnostic; correctness first")
println("CPU: ", Sys.CPU_NAME)
run_diagnostics()
run_timing()
