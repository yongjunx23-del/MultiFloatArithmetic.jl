using MultiFloatArithmetic
using MultiFloats
using Random

const E = MultiFloatArithmetic.Experimental
const T = MultiFloat{Float64,5}

function pow2_fma5_bench(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

function dense_fma5_bench(; emin=-70, emax=70)
    numerator = BigInt(0)
    for _ in 1:5
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64 * 5 - 1)
    q = numerator // (BigInt(1) << (64 * 5))
    q *= pow2_fma5_bench(rand(emin:emax))
    q = rand(Bool) ? q : -q
    return setprecision(BigFloat, 1536) do
        T(BigFloat(q))
    end
end

mutable struct Fma5Stats
    n::Int
    oracle_mismatches::Int
    normalization_failures::Int
    symmetry_failures::Int
    composition_oracle_mismatches::Int
    nonzero_tail::Int
    max_operand_constant::BigFloat
    max_result_relative::BigFloat
end
Fma5Stats() = Fma5Stats(0, 0, 0, 0, 0, 0, BigFloat(0), BigFloat(0))

function evaluate!(stats::Fma5Stats, x::T, y::T, c::T)
    direct = E.fma5_safe(x, y, c)
    oracle = E.reference_fma(x, y, c)
    composed = E.add_safe(E.mul_safe(x, y), c)

    qx = Rational{BigInt}(x)
    qy = Rational{BigInt}(y)
    qc = Rational{BigInt}(c)
    qz = Rational{BigInt}(direct)
    qexact = qx * qy + qc
    err = abs(qz - qexact)

    stats.n += 1
    stats.oracle_mismatches += direct !== oracle
    stats.normalization_failures += !MultiFloats.isnormalized(direct)
    stats.symmetry_failures += direct !== E.fma5_safe(y, x, c)
    stats.composition_oracle_mismatches += composed !== oracle
    stats.nonzero_tail += !iszero(err)

    if !iszero(err)
        setprecision(BigFloat, 768) do
            xy = BigFloat(qx * qy)
            cb = BigFloat(qc)
            scale = abs(xy) + abs(cb)
            if !iszero(scale)
                constant = BigFloat(err) / (BigFloat(2)^(-53 * 5) * scale)
                stats.max_operand_constant = max(stats.max_operand_constant, constant)
            end
            result_scale = abs(BigFloat(qexact))
            if !iszero(result_scale)
                stats.max_result_relative = max(
                    stats.max_result_relative,
                    BigFloat(err) / result_scale,
                )
            end
        end
    end
    return nothing
end

function report(label, stats::Fma5Stats)
    println(label)
    println("  cases=$(stats.n)")
    println("  direct_oracle_mismatches=$(stats.oracle_mismatches)")
    println("  normalization_failures=$(stats.normalization_failures)")
    println("  x/y_symmetry_failures=$(stats.symmetry_failures)")
    println("  composed_mul_then_add_oracle_mismatches=$(stats.composition_oracle_mismatches)")
    println("  nonzero_discarded_tail=$(stats.nonzero_tail)")
    println("  max |err|/(u^5*(|xy|+|c|)) = $(stats.max_operand_constant)")
    println("  max result-relative error = $(stats.max_result_relative)")
end

function assert_gate(stats::Fma5Stats)
    @assert stats.oracle_mismatches == 0
    @assert stats.normalization_failures == 0
    @assert stats.symmetry_failures == 0
    @assert stats.max_operand_constant <= 1
end

function cancellation_case(depth)
    x = dense_fma5_bench(emin=-50, emax=50)
    y = dense_fma5_bench(emin=-50, emax=50)
    qxy = Rational{BigInt}(x) * Rational{BigInt}(y)
    lead_exp = exponent(x._limbs[1]) + exponent(y._limbs[1])
    delta = pow2_fma5_bench(lead_exp - depth)
    c = setprecision(BigFloat, 1536) do
        T(BigFloat(-qxy + (rand(Bool) ? delta : -delta)))
    end
    return x, y, c
end

function min_time(f; samples=3)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function run_diagnostics()
    Random.seed!(0xf5a5_2026)

    ordinary = Fma5Stats()
    for _ in 1:200
        evaluate!(ordinary, dense_fma5_bench(), dense_fma5_bench(), dense_fma5_bench())
    end
    report("dense ordinary", ordinary)
    assert_gate(ordinary)

    scaled = Fma5Stats()
    for _ in 1:200
        evaluate!(scaled,
                  dense_fma5_bench(emin=-100, emax=100),
                  dense_fma5_bench(emin=-100, emax=100),
                  dense_fma5_bench(emin=-100, emax=100))
    end
    report("scaled", scaled)
    assert_gate(scaled)

    cancellation = Fma5Stats()
    for depth in (80, 160, 230, 270, 320)
        for _ in 1:30
            x, y, c = cancellation_case(depth)
            evaluate!(cancellation, x, y, c)
        end
    end
    report("destructive cancellation", cancellation)
    assert_gate(cancellation)
end

function run_timing(; n=120)
    Random.seed!(0xf5a0_2026)
    xs = [dense_fma5_bench(emin=-40, emax=40) for _ in 1:n]
    ys = [dense_fma5_bench(emin=-40, emax=40) for _ in 1:n]
    cs = [dense_fma5_bench(emin=-40, emax=40) for _ in 1:n]
    out = Vector{T}(undef, n)

    direct!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.fma5_safe(xs[i], ys[i], cs[i])
        end
        out
    end
    composed!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.add_safe(E.mul_safe(xs[i], ys[i]), cs[i])
        end
        out
    end

    direct!(); composed!()
    td = min_time(direct!)
    tc = min_time(composed!)
    println("Float64x5 direct safe FMA: $(round(td*1e3; digits=3)) ms for $(n) cases")
    println("Float64x5 safe mul-then-add: $(round(tc*1e3; digits=3)) ms for $(n) cases")
    println("composition/direct time ratio = $(round(tc/td; digits=3))x")
end

println("M5 Float64x5 direct-FMA diagnostic; correctness first")
println("CPU: ", Sys.CPU_NAME)
run_diagnostics()
run_timing()
