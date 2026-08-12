using MultiFloatArithmetic
using MultiFloats
using Random

const E = MultiFloatArithmetic.Experimental

function wide_randN(::Type{T}; emin=-350, emax=350) where {T}
    x = rand(Bool) ? rand(T) : -rand(T)
    e = rand(emin:emax)
    return T(ldexp(big(x), e))
end

function pow2q(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

mutable struct AddStats
    n::Int
    mismatches::Int
    normalization_failures::Int
    commutativity_failures::Int
    nonzero_tail::Int
    max_operand_constant::BigFloat
end
AddStats() = AddStats(0, 0, 0, 0, 0, BigFloat(0))

function evaluate!(stats::AddStats, x::MultiFloat{Float64,N}, y::MultiFloat{Float64,N}) where {N}
    z = E.add_safe(x, y)
    oracle = E.reference_add(x, y)
    qx = Rational{BigInt}(x)
    qy = Rational{BigInt}(y)
    qz = Rational{BigInt}(z)
    qexact = qx + qy
    err = abs(qz - qexact)

    stats.n += 1
    stats.mismatches += z !== oracle
    stats.normalization_failures += !MultiFloats.isnormalized(z)
    stats.commutativity_failures += z !== E.add_safe(y, x)
    stats.nonzero_tail += !iszero(err)

    if !iszero(err)
        setprecision(BigFloat, 768) do
            scale = abs(BigFloat(qx)) + abs(BigFloat(qy))
            if !iszero(scale)
                constant = BigFloat(err) / (BigFloat(2)^(-53N) * scale)
                stats.max_operand_constant = max(stats.max_operand_constant, constant)
            end
        end
    end
    return nothing
end

function assert_gate(stats::AddStats)
    @assert stats.mismatches == 0
    @assert stats.normalization_failures == 0
    @assert stats.commutativity_failures == 0
    # Width-specific empirical regression gate. C=1 is intentionally loose and
    # remains non-theorem evidence until M3 formal tail analysis is complete.
    @assert stats.max_operand_constant <= 1
end

function report(N, label, stats)
    println("Float64x$(N) $(label)")
    println("  cases=$(stats.n), oracle_mismatches=$(stats.mismatches), ",
            "normalization_failures=$(stats.normalization_failures), ",
            "commutativity_failures=$(stats.commutativity_failures)")
    println("  nonzero_discarded_tail=$(stats.nonzero_tail)")
    println("  max |err|/(u^$(N)*(|x|+|y|)) = $(stats.max_operand_constant)")
end

function run_width(N)
    T = MultiFloat{Float64,N}
    Random.seed!(0xad00_2026 + N)

    ordinary = AddStats()
    for _ in 1:600
        evaluate!(ordinary,
                  rand(Bool) ? rand(T) : -rand(T),
                  rand(Bool) ? rand(T) : -rand(T))
    end
    report(N, "ordinary", ordinary)
    assert_gate(ordinary)

    wide = AddStats()
    for _ in 1:600
        evaluate!(wide, wide_randN(T), wide_randN(T))
    end
    report(N, "wide exponent", wide)
    assert_gate(wide)

    cancellation = AddStats()
    for depth in (100, 53N - 60, 53N - 10, 53N + 30)
        for _ in 1:50
            x = wide_randN(T; emin=-150, emax=150)
            qx = Rational{BigInt}(x)
            lead = x._limbs[1]
            e = iszero(lead) ? 0 : exponent(lead)
            delta = pow2q(e - depth)
            y = T(-qx + (rand(Bool) ? delta : -delta))
            evaluate!(cancellation, x, y)
        end
    end
    report(N, "near cancellation", cancellation)
    assert_gate(cancellation)
end

println("M3 Float64x6-x8 safe-addition diagnostics")
println("CPU: ", Sys.CPU_NAME)
for N in 6:8
    run_width(N)
end
