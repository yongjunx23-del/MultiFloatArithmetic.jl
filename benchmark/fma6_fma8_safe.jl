using MultiFloatArithmetic
using MultiFloats
using Random

const E = MultiFloatArithmetic.Experimental

function pow2_fmaN_bench(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

function dense_fmaN_bench(::Type{T}, N; emin=-50, emax=50) where {T}
    numerator = BigInt(0)
    for _ in 1:N
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64N - 1)
    q = numerator // (BigInt(1) << 64N)
    q *= pow2_fmaN_bench(rand(emin:emax))
    q = rand(Bool) ? q : -q
    return setprecision(BigFloat, 2048) do
        T(BigFloat(q))
    end
end

mutable struct FmaNStats
    n::Int
    direct_oracle_mismatches::Int
    normalization_failures::Int
    symmetry_failures::Int
    composition_oracle_mismatches::Int
    nonzero_tail::Int
    max_operand_constant::BigFloat
    max_result_relative::BigFloat
end
FmaNStats() = FmaNStats(0, 0, 0, 0, 0, 0, BigFloat(0), BigFloat(0))

function evaluate!(
    stats::FmaNStats,
    x::MultiFloat{Float64,N},
    y::MultiFloat{Float64,N},
    c::MultiFloat{Float64,N},
) where {N}
    direct = E.fma_safe(x, y, c)
    oracle = E.reference_fma(x, y, c)
    composed = E.add_safe(E.mul_safe(x, y), c)

    qx = Rational{BigInt}(x)
    qy = Rational{BigInt}(y)
    qc = Rational{BigInt}(c)
    qz = Rational{BigInt}(direct)
    exact = qx * qy + qc
    err = abs(qz - exact)

    stats.n += 1
    stats.direct_oracle_mismatches += direct !== oracle
    stats.normalization_failures += !MultiFloats.isnormalized(direct)
    stats.symmetry_failures += direct !== E.fma_safe(y, x, c)
    stats.composition_oracle_mismatches += composed !== oracle
    stats.nonzero_tail += !iszero(err)

    if !iszero(err)
        setprecision(BigFloat, 1152) do
            xy = BigFloat(qx * qy)
            cb = BigFloat(qc)
            scale = abs(xy) + abs(cb)
            if !iszero(scale)
                constant = BigFloat(err) / (BigFloat(2)^(-53N) * scale)
                stats.max_operand_constant = max(stats.max_operand_constant, constant)
            end
            result_scale = abs(BigFloat(exact))
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

function report(N, label, stats::FmaNStats)
    println("Float64x$(N) $(label)")
    println("  cases=$(stats.n), direct_oracle_mismatches=$(stats.direct_oracle_mismatches), ",
            "normalization_failures=$(stats.normalization_failures), ",
            "x/y_symmetry_failures=$(stats.symmetry_failures)")
    println("  composed_mul_then_add_oracle_mismatches=$(stats.composition_oracle_mismatches)")
    println("  nonzero_discarded_tail=$(stats.nonzero_tail)")
    println("  max |err|/(u^$(N)*(|xy|+|c|)) = $(stats.max_operand_constant)")
    println("  max result-relative error = $(stats.max_result_relative)")
end

function assert_gate(stats::FmaNStats)
    @assert stats.direct_oracle_mismatches == 0
    @assert stats.normalization_failures == 0
    @assert stats.symmetry_failures == 0
    @assert stats.max_operand_constant <= 1
end

function cancellation_case(::Type{T}, N, depth) where {T}
    x = dense_fmaN_bench(T, N; emin=-40, emax=40)
    y = dense_fmaN_bench(T, N; emin=-40, emax=40)
    qxy = Rational{BigInt}(x) * Rational{BigInt}(y)
    lead_exp = exponent(x._limbs[1]) + exponent(y._limbs[1])
    delta = pow2_fmaN_bench(lead_exp - depth)
    c = setprecision(BigFloat, 2048) do
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

function run_width(N)
    T = MultiFloat{Float64,N}
    Random.seed!(0xf600_2026 + N)
    ncase = N == 6 ? 50 : N == 7 ? 35 : 25

    ordinary = FmaNStats()
    for _ in 1:ncase
        evaluate!(ordinary,
                  dense_fmaN_bench(T, N),
                  dense_fmaN_bench(T, N),
                  dense_fmaN_bench(T, N))
    end
    report(N, "dense ordinary", ordinary)
    assert_gate(ordinary)

    scaled = FmaNStats()
    for _ in 1:ncase
        evaluate!(scaled,
                  dense_fmaN_bench(T, N; emin=-80, emax=80),
                  dense_fmaN_bench(T, N; emin=-80, emax=80),
                  dense_fmaN_bench(T, N; emin=-80, emax=80))
    end
    report(N, "scaled", scaled)
    assert_gate(scaled)

    cancellation = FmaNStats()
    ncancel = N == 6 ? 12 : N == 7 ? 8 : 6
    for depth in (100, 53N - 80, 53N - 20, 53N + 30), _ in 1:ncancel
        x, y, c = cancellation_case(T, N, depth)
        evaluate!(cancellation, x, y, c)
    end
    report(N, "destructive cancellation", cancellation)
    assert_gate(cancellation)

    ntiming = N == 6 ? 20 : N == 7 ? 12 : 6
    xs = [dense_fmaN_bench(T, N; emin=-35, emax=35) for _ in 1:ntiming]
    ys = [dense_fmaN_bench(T, N; emin=-35, emax=35) for _ in 1:ntiming]
    cs = [dense_fmaN_bench(T, N; emin=-35, emax=35) for _ in 1:ntiming]
    out = Vector{T}(undef, ntiming)

    direct!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = E.fma_safe(xs[i], ys[i], cs[i])
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
    println("Float64x$(N) direct safe FMA: $(round(td*1e3; digits=3)) ms for $(ntiming) cases")
    println("Float64x$(N) safe mul-then-add: $(round(tc*1e3; digits=3)) ms for $(ntiming) cases")
    println("Float64x$(N) composition/direct time ratio = $(round(tc/td; digits=3))x")
end

println("M5 Float64x6-x8 direct-FMA diagnostics")
println("CPU: ", Sys.CPU_NAME)
for N in 6:8
    run_width(N)
end
