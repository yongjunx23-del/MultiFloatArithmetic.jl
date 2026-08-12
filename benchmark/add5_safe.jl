using MultiFloatArithmetic
using MultiFloats
using Random

const E = MultiFloatArithmetic.Experimental
const T = MultiFloat{Float64,5}

@inline exact(x::T) = Rational{BigInt}(x)

function wide_rand5(; emin=-400, emax=400)
    x = rand(Bool) ? rand(T) : -rand(T)
    e = rand(emin:emax)
    return T(ldexp(big(x), e))
end

function cancellation_pair(bits::Int)
    x = wide_rand5(emin=-200, emax=200)
    qx = exact(x)
    lead = x._limbs[1]
    e = iszero(lead) ? 0 : exponent(lead)
    delta = if e >= bits
        (BigInt(1) << (e - bits)) // BigInt(1)
    else
        BigInt(1) // (BigInt(1) << (bits - e))
    end
    y = T(-qx + (rand(Bool) ? delta : -delta))
    return x, y
end

function min_time(f; samples=7)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

mutable struct Stats
    n::Int
    mismatches::Int
    normalization_failures::Int
    commutativity_failures::Int
    nonzero_tail::Int
    max_operand_constant::BigFloat
    max_result_relative::BigFloat
end
Stats() = Stats(0, 0, 0, 0, 0, BigFloat(0), BigFloat(0))

function evaluate!(stats::Stats, x::T, y::T)
    z = E.add5_safe(x, y)
    oracle = E.reference_add(x, y)
    qx = exact(x)
    qy = exact(y)
    qz = exact(z)
    qexact = qx + qy
    err = abs(qz - qexact)

    stats.n += 1
    stats.mismatches += z !== oracle
    stats.normalization_failures += !MultiFloats.isnormalized(z)
    stats.commutativity_failures += z !== E.add5_safe(y, x)
    stats.nonzero_tail += !iszero(err)

    if !iszero(err)
        setprecision(BigFloat, 512) do
            err_b = BigFloat(err)
            scale = abs(BigFloat(qx)) + abs(BigFloat(qy))
            u5 = BigFloat(2)^(-53 * 5)
            if !iszero(scale)
                stats.max_operand_constant = max(
                    stats.max_operand_constant,
                    err_b / (u5 * scale),
                )
            end
            result_scale = abs(BigFloat(qexact))
            if !iszero(result_scale)
                stats.max_result_relative = max(
                    stats.max_result_relative,
                    err_b / result_scale,
                )
            end
        end
    end
    return nothing
end

function report(label, stats::Stats)
    println(label)
    println("  cases=$(stats.n)")
    println("  oracle_mismatches=$(stats.mismatches)")
    println("  normalization_failures=$(stats.normalization_failures)")
    println("  commutativity_failures=$(stats.commutativity_failures)")
    println("  nonzero_discarded_tail=$(stats.nonzero_tail)")
    println("  max |err|/(u^5*(|x|+|y|)) = ", stats.max_operand_constant)
    println("  max result-relative error = ", stats.max_result_relative)
end

function assert_acceptance_gate(stats::Stats)
    @assert stats.mismatches == 0
    @assert stats.normalization_failures == 0
    @assert stats.commutativity_failures == 0
    # Empirical regression gate only, not a proved theorem. The first accepted
    # corpus measured about 0.05284, so C=1 leaves nearly 19x headroom while
    # still catching an order-of-magnitude tail regression.
    @assert stats.max_operand_constant <= 1
end

function run_diagnostics()
    Random.seed!(0xad05_2026)

    ordinary = Stats()
    for _ in 1:1_000
        x = rand(Bool) ? rand(T) : -rand(T)
        y = rand(Bool) ? rand(T) : -rand(T)
        evaluate!(ordinary, x, y)
    end
    report("ordinary", ordinary)

    wide = Stats()
    for _ in 1:1_000
        evaluate!(wide, wide_rand5(), wide_rand5())
    end
    report("wide exponent", wide)

    cancellation = Stats()
    for bits in (100, 180, 240, 280)
        for _ in 1:125
            x, y = cancellation_pair(bits)
            evaluate!(cancellation, x, y)
        end
    end
    report("near cancellation", cancellation)

    assert_acceptance_gate(ordinary)
    assert_acceptance_gate(wide)
    assert_acceptance_gate(cancellation)
end

function run_timing(; n=5_000)
    Random.seed!(0xadd5_2026)
    xs = [rand(Bool) ? rand(T) : -rand(T) for _ in 1:n]
    ys = [rand(Bool) ? rand(T) : -rand(T) for _ in 1:n]
    out = Vector{T}(undef, n)

    safe!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.add5_safe(xs[i], ys[i])
        end
        out
    end

    safe!()
    ts = min_time(safe!)
    println("Float64x5 safe addition: $(round(ts*1e3; digits=3)) ms for $(n) cases")
end

println("M3 Float64x5 safe-addition diagnostic; correctness first")
println("CPU: ", Sys.CPU_NAME)
run_diagnostics()
run_timing()
