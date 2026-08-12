using MultiFloatArithmetic
using MultiFloats
using Random

const E = MultiFloatArithmetic.Experimental

function pow2_mulN_bench(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

function dense_mulN_bench(::Type{T}, N; emin=-50, emax=50) where {T}
    numerator = BigInt(0)
    for _ in 1:N
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64N - 1)
    q = numerator // (BigInt(1) << 64N)
    q *= pow2_mulN_bench(rand(emin:emax))
    q = rand(Bool) ? q : -q
    return setprecision(BigFloat, 2048) do
        T(BigFloat(q))
    end
end

mutable struct MulNStats
    n::Int
    mismatches::Int
    normalization_failures::Int
    commutativity_failures::Int
    nonzero_tail::Int
    max_relative_constant::BigFloat
end
MulNStats() = MulNStats(0, 0, 0, 0, 0, BigFloat(0))

function evaluate!(stats::MulNStats, x::MultiFloat{Float64,N}, y::MultiFloat{Float64,N}) where {N}
    z = E.mul_safe(x, y)
    oracle = E.reference_mul(x, y)
    qx = Rational{BigInt}(x)
    qy = Rational{BigInt}(y)
    qz = Rational{BigInt}(z)
    exact = qx * qy
    err = abs(qz - exact)

    stats.n += 1
    stats.mismatches += z !== oracle
    stats.normalization_failures += !MultiFloats.isnormalized(z)
    stats.commutativity_failures += z !== E.mul_safe(y, x)
    stats.nonzero_tail += !iszero(err)

    if !iszero(err) && !iszero(exact)
        setprecision(BigFloat, 1024) do
            constant = BigFloat(err) / (BigFloat(2)^(-53N) * abs(BigFloat(exact)))
            stats.max_relative_constant = max(stats.max_relative_constant, constant)
        end
    end
    return nothing
end

function report(N, label, stats::MulNStats)
    println("Float64x$(N) $(label)")
    println("  cases=$(stats.n), oracle_mismatches=$(stats.mismatches), ",
            "normalization_failures=$(stats.normalization_failures), ",
            "commutativity_failures=$(stats.commutativity_failures)")
    println("  nonzero_discarded_tail=$(stats.nonzero_tail)")
    println("  max |err|/(u^$(N)*|x*y|) = $(stats.max_relative_constant)")
end

function assert_gate(stats::MulNStats)
    @assert stats.mismatches == 0
    @assert stats.normalization_failures == 0
    @assert stats.commutativity_failures == 0
    @assert stats.max_relative_constant <= 1
end

function min_time(f; samples=3)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function run_width(N)
    T = MultiFloat{Float64,N}
    Random.seed!(0x4d00_2026 + N)
    ncase = N == 6 ? 80 : N == 7 ? 60 : 40

    ordinary = MulNStats()
    for _ in 1:ncase
        evaluate!(ordinary, dense_mulN_bench(T, N), dense_mulN_bench(T, N))
    end
    report(N, "dense ordinary", ordinary)
    assert_gate(ordinary)

    scaled = MulNStats()
    for _ in 1:ncase
        evaluate!(scaled,
                  dense_mulN_bench(T, N; emin=-90, emax=90),
                  dense_mulN_bench(T, N; emin=-90, emax=90))
    end
    report(N, "scaled", scaled)
    assert_gate(scaled)

    boundary = MulNStats()
    for ex in (-80, 0, 80), ey in (-80, 0, 80)
        x = T(BigFloat(pow2_mulN_bench(ex)))
        y = T(BigFloat(pow2_mulN_bench(ey)))
        evaluate!(boundary, x, y)
    end
    report(N, "power-of-two boundaries", boundary)
    assert_gate(boundary)

    ntiming = N == 6 ? 40 : N == 7 ? 25 : 12
    xs = [dense_mulN_bench(T, N; emin=-40, emax=40) for _ in 1:ntiming]
    ys = [dense_mulN_bench(T, N; emin=-40, emax=40) for _ in 1:ntiming]
    out = Vector{T}(undef, ntiming)
    work!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.mul_safe(xs[i], ys[i])
        end
        out
    end
    work!()
    t = min_time(work!)
    println("Float64x$(N) safe multiplication: $(round(t*1e3; digits=3)) ms for $(ntiming) cases")
end

println("M4 Float64x6-x8 safe-multiplication diagnostics")
println("CPU: ", Sys.CPU_NAME)
for N in 6:8
    run_width(N)
end
