using MultiFloatArithmetic
using MultiFloats
using Random
using Printf

const E = MultiFloatArithmetic.Experimental

# Requested comparison pairs. A Float64 limb contributes at most 53 significand
# bits, so these are deliberately NOT equal-precision pairs:
#   Float64x2 ~= 106 bits vs BigFloat128
#   Float64x4 ~= 212 bits vs BigFloat256
#   Float64x8 ~= 424 bits vs BigFloat512
const PAIRS = ((2, 128), (4, 256), (8, 512))

function pow2_cmp(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) :
        BigInt(1) // (BigInt(1) << (-e))
end

function source_value(; words=8, emin=-20, emax=20)
    numerator = BigInt(0)
    for _ in 1:words
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64words - 1)
    q = numerator // (BigInt(1) << 64words)
    return q * pow2_cmp(rand(emin:emax))
end

function make_sources(n)
    xs = [source_value() for _ in 1:n]
    ys = [source_value() for _ in 1:n]
    cs = [source_value() for _ in 1:n]
    return xs, ys, cs
end

function make_bigfloats(qs, bits)
    return setprecision(BigFloat, bits) do
        [BigFloat(q) for q in qs]
    end
end

function make_multifloats(qs, N)
    T = MultiFloat{Float64,N}
    return setprecision(BigFloat, 2048) do
        [T(BigFloat(q)) for q in qs]
    end
end

function relative_error(qhat::Rational{BigInt}, qref::Rational{BigInt})
    iszero(qref) && return iszero(qhat) ? BigFloat(0) : BigFloat(Inf)
    return setprecision(BigFloat, 2048) do
        BigFloat(abs(qhat - qref)) / abs(BigFloat(qref))
    end
end

function effective_bits(err::BigFloat)
    iszero(err) && return Inf
    return Float64(-log2(err))
end

function max_errors_bigfloat(xs, ys, cs, qx, qy, qc, bits)
    max_add = BigFloat(0)
    max_mul = BigFloat(0)
    max_ma = BigFloat(0)
    setprecision(BigFloat, bits) do
        for i in eachindex(xs)
            za = xs[i] + ys[i]
            zm = xs[i] * ys[i]
            zma = xs[i] * ys[i] + cs[i]
            max_add = max(max_add,
                relative_error(Rational{BigInt}(za), qx[i] + qy[i]))
            max_mul = max(max_mul,
                relative_error(Rational{BigInt}(zm), qx[i] * qy[i]))
            max_ma = max(max_ma,
                relative_error(Rational{BigInt}(zma), qx[i] * qy[i] + qc[i]))
        end
    end
    return max_add, max_mul, max_ma
end

function max_errors_multifloat(xs, ys, cs, qx, qy, qc, N)
    max_add = BigFloat(0)
    max_mul = BigFloat(0)
    max_ma_original = N <= 4 ? BigFloat(0) : nothing
    max_ma_project = N <= 4 ? BigFloat(0) : nothing

    for i in eachindex(xs)
        if N <= 4
            za = xs[i] + ys[i]               # original MultiFloats.jl
            zm = xs[i] * ys[i]               # original MultiFloats.jl
            zma = xs[i] * ys[i] + cs[i]      # original MultiFloats.jl
            zf = fma_fast(xs[i], ys[i], cs[i])
            max_ma_original = max(max_ma_original,
                relative_error(Rational{BigInt}(zma), qx[i] * qy[i] + qc[i]))
            max_ma_project = max(max_ma_project,
                relative_error(Rational{BigInt}(zf), qx[i] * qy[i] + qc[i]))
        else
            za = E.add_safe(xs[i], ys[i])
            zm = E.mul_safe(xs[i], ys[i])
        end

        max_add = max(max_add,
            relative_error(Rational{BigInt}(za), qx[i] + qy[i]))
        max_mul = max(max_mul,
            relative_error(Rational{BigInt}(zm), qx[i] * qy[i]))
    end
    return max_add, max_mul, max_ma_original, max_ma_project
end

function best_ns_per_op(work!, nops; samples=3)
    work!()
    best = Inf
    for _ in 1:samples
        GC.gc()
        t = @elapsed work!()
        best = min(best, t)
    end
    return best * 1e9 / nops
end

function time_bigfloat(xs, ys, cs, bits, op; reps)
    out = Vector{BigFloat}(undef, length(xs))
    work! = if op === :add
        () -> setprecision(BigFloat, bits) do
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = xs[i] + ys[i]
            end
        end
    elseif op === :mul
        () -> setprecision(BigFloat, bits) do
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = xs[i] * ys[i]
            end
        end
    else
        () -> setprecision(BigFloat, bits) do
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = xs[i] * ys[i] + cs[i]
            end
        end
    end
    return best_ns_per_op(work!, reps * length(xs))
end

function time_multifloat(xs, ys, cs, N, op; reps)
    T = MultiFloat{Float64,N}
    out = Vector{T}(undef, length(xs))
    work! = if N <= 4 && op === :add
        () -> begin
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = xs[i] + ys[i]
            end
        end
    elseif N <= 4 && op === :mul
        () -> begin
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = xs[i] * ys[i]
            end
        end
    elseif N <= 4 && op === :muladd_original
        () -> begin
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = xs[i] * ys[i] + cs[i]
            end
        end
    elseif N <= 4 && op === :fma_fast
        () -> begin
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = fma_fast(xs[i], ys[i], cs[i])
            end
        end
    elseif N == 8 && op === :add
        () -> begin
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = E.add8_safe(xs[i], ys[i])
            end
        end
    elseif N == 8 && op === :mul
        () -> begin
            for _ in 1:reps, i in eachindex(xs)
                @inbounds out[i] = E.mul8_safe(xs[i], ys[i])
            end
        end
    else
        error("unsupported timing combination N=$(N), op=$(op)")
    end
    return best_ns_per_op(work!, reps * length(xs))
end

function print_error(label, err)
    @printf("    %-24s max rel = %.4e, observed bits ~= %.1f\n",
        label, Float64(err), effective_bits(err))
end

function run_pair(N, bits, qx, qy, qc; precision_n=96, timing_n=256)
    qxp = qx[1:precision_n]
    qyp = qy[1:precision_n]
    qcp = qc[1:precision_n]

    bx = make_bigfloats(qxp, bits)
    by = make_bigfloats(qyp, bits)
    bc = make_bigfloats(qcp, bits)
    mx = make_multifloats(qxp, N)
    my = make_multifloats(qyp, N)
    mc = make_multifloats(qcp, N)

    badd, bmul, bma = max_errors_bigfloat(bx, by, bc, qxp, qyp, qcp, bits)
    madd, mmul, mma_original, mma_project =
        max_errors_multifloat(mx, my, mc, qxp, qyp, qcp, N)

    println()
    println("BigFloat$(bits) vs Float64x$(N)")
    println("  nominal significand capacity: BigFloat=$(bits) bits; Float64x$(N) ~= $(53N) bits")
    println("  precision: end-to-end error against shared exact rational source values")
    print_error("BigFloat add", badd)
    print_error("MultiFloat add", madd)
    print_error("BigFloat mul", bmul)
    print_error("MultiFloat mul", mmul)
    if N <= 4
        print_error("BigFloat mul+add", bma)
        print_error("original MF mul+add", mma_original)
        print_error("project fma_fast", mma_project)
    end

    qxt = qx[1:timing_n]
    qyt = qy[1:timing_n]
    qct = qc[1:timing_n]
    bxt = make_bigfloats(qxt, bits)
    byt = make_bigfloats(qyt, bits)
    bct = make_bigfloats(qct, bits)
    mxt = make_multifloats(qxt, N)
    myt = make_multifloats(qyt, N)
    mct = make_multifloats(qct, N)

    breps = bits == 128 ? 100 : bits == 256 ? 60 : 25
    mreps = N == 2 ? 200 : N == 4 ? 100 : 8
    mmulreps = N == 8 ? 2 : mreps

    badd_ns = time_bigfloat(bxt, byt, bct, bits, :add; reps=breps)
    bmul_ns = time_bigfloat(bxt, byt, bct, bits, :mul; reps=breps)
    madd_ns = time_multifloat(mxt, myt, mct, N, :add; reps=mreps)
    mmul_ns = time_multifloat(mxt, myt, mct, N, :mul; reps=mmulreps)

    println("  speed: best hosted-runner scalar time, ns/op; lower is faster")
    @printf("    add: BigFloat=%9.1f  MultiFloat=%9.1f  BigFloat/MultiFloat=%6.2fx\n",
        badd_ns, madd_ns, badd_ns / madd_ns)
    @printf("    mul: BigFloat=%9.1f  MultiFloat=%9.1f  BigFloat/MultiFloat=%6.2fx\n",
        bmul_ns, mmul_ns, bmul_ns / mmul_ns)

    if N <= 4
        bma_ns = time_bigfloat(bxt, byt, bct, bits, :muladd; reps=breps)
        oma_ns = time_multifloat(mxt, myt, mct, N, :muladd_original; reps=mreps)
        fma_ns = time_multifloat(mxt, myt, mct, N, :fma_fast; reps=mreps)
        @printf("    mul+add: BigFloat=%9.1f  original MultiFloat=%9.1f  project fma_fast=%9.1f ns/op\n",
            bma_ns, oma_ns, fma_ns)
    else
        println("    original MultiFloats.jl x8 arithmetic: unavailable upstream; x8 rows use this project's safe baselines")
    end
end

println("BigFloat versus MultiFloat scalar precision/speed comparison")
println("CPU: ", Sys.CPU_NAME)
println("MultiFloats version: 3.2.6 (pinned by project manifest)")
println("Precision note: requested pairs give BigFloat about 20.8% more nominal significand bits than Float64xN.")

Random.seed!(0xb16f_2026)
qx, qy, qc = make_sources(256)
for (N, bits) in PAIRS
    run_pair(N, bits, qx, qy, qc)
end
