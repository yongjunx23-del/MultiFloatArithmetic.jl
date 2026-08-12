using MultiFloatArithmetic
using MultiFloats
using Random
using Printf

include(joinpath(@__DIR__, "..", "audit", "paper2607_v4.jl"))
using .Paper2607Audit

const T2P = Float64x2
const T4P = Float64x4

function minimum_time(f; samples=7)
    f()
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

paper_exact(x) = Rational{BigInt}(x)

function signed_rand(::Type{T}) where {T}
    x = rand(T)
    return MultiFloats.renormalize(rand(Bool) ? x : -x)
end

function wide_rand(::Type{T}; emin=-350, emax=350) where {T}
    return MultiFloats.renormalize(ldexp(signed_rand(T), rand(emin:emax)))
end

function pack_vectors(::Val{W}, scalars) where {W}
    V = MultiFloatVec{W,Float64,4}
    nvec = length(scalars) ÷ W
    out = Vector{V}(undef, nvec)
    @inbounds for i in 1:nvec
        base = W * (i - 1)
        out[i] = V(ntuple(lane -> scalars[base + lane], Val(W)))
    end
    return out
end

function audit_dw(; n=20_000)
    Random.seed!(0x2607_d002)
    xs = [signed_rand(T2P) for _ in 1:n]
    ys = [signed_rand(T2P) for _ in 1:n]
    cs = [signed_rand(T2P) for _ in 1:n]
    paper = similar(xs)
    current = similar(xs)
    upstream = similar(xs)

    paper!() = begin
        @inbounds for i in eachindex(xs)
            paper[i] = paper_dw_v4(xs[i], ys[i], cs[i])
        end
        paper
    end
    current!() = begin
        @inbounds for i in eachindex(xs)
            current[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        current
    end
    upstream!() = begin
        @inbounds for i in eachindex(xs)
            upstream[i] = xs[i] * ys[i] + cs[i]
        end
        upstream
    end

    paper!(); current!(); upstream!()
    mismatches = count(i -> paper[i] !== current[i], eachindex(paper))
    tp = minimum_time(paper!)
    tc = minimum_time(current!)
    tu = minimum_time(upstream!)

    println("DW / Float64x2 — arXiv v4 Algorithm 1")
    println("  paper-v4 vs project fma_fast bitwise mismatches = $mismatches/$n")
    @printf("  paper-v4 17-flop FMA = %.3f ms / %d\n", 1e3*tp, n)
    @printf("  current fma_fast     = %.3f ms / %d\n", 1e3*tc, n)
    @printf("  upstream x*y+c       = %.3f ms / %d\n", 1e3*tu, n)
    @printf("  upstream/paper time ratio = %.3fx\n", tu/tp)
end

function audit_qw_corpus(label, xs, ys, cs; exact_sample=1_000)
    n = length(xs)
    v1_bad = 0
    v4_bad = 0
    safe_bad = 0
    f2s_violations = zeros(Int, 10)
    v4_safe_bitwise = 0
    exact_v1_v4_mismatch = 0
    exact_v4_safe_mismatch = 0
    direct_better = 0
    upstream_better = 0
    equal_error = 0

    setprecision(BigFloat, 1024) do
        @inbounds for i in eachindex(xs)
            x, y, c = xs[i], ys[i], cs[i]
            v1 = paper_qw_v1_2pass(x, y, c)
            v4 = paper_qw_v4_5pass(x, y, c)
            safe = fma_fast(x, y, c)
            upstream = x * y + c

            v1_bad += !MultiFloats.isnormalized(v1)
            v4_bad += !MultiFloats.isnormalized(v4)
            safe_bad += !MultiFloats.isnormalized(safe)
            checks = paper_qw_v4_fast_preconditions(x, y, c)
            for j in eachindex(checks)
                f2s_violations[j] += !checks[j]
            end
            v4_safe_bitwise += (v4 !== safe)

            if i <= exact_sample
                rv1 = paper_exact(v1)
                rv4 = paper_exact(v4)
                rsafe = paper_exact(safe)
                exact_v1_v4_mismatch += (rv1 != rv4)
                exact_v4_safe_mismatch += (rv4 != rsafe)

                truth = paper_exact(x) * paper_exact(y) + paper_exact(c)
                ev4 = abs(rv4 - truth)
                eup = abs(paper_exact(upstream) - truth)
                if ev4 < eup
                    direct_better += 1
                elseif eup < ev4
                    upstream_better += 1
                else
                    equal_error += 1
                end
            end
        end
    end

    println(label)
    println("  cases=$n")
    println("  non-normalized: v1-2pass=$v1_bad, v4-5pass=$v4_bad, current-safe=$safe_bad")
    println("  v4-vs-current bitwise mismatches=$v4_safe_bitwise/$n")
    println("  paper FastTwoSum exponent-contract violations by gate = $f2s_violations")
    println("  exact component-sum mismatches (first $exact_sample): v1-v4=$exact_v1_v4_mismatch, v4-safe=$exact_v4_safe_mismatch")
    println("  exact-error comparison vs upstream (first $exact_sample): direct_better=$direct_better, upstream_better=$upstream_better, equal=$equal_error")
end

function audit_qw_scalar_speed(; n=20_000)
    Random.seed!(0x2607_4004)
    xs = [signed_rand(T4P) for _ in 1:n]
    ys = [signed_rand(T4P) for _ in 1:n]
    cs = [signed_rand(T4P) for _ in 1:n]
    v1 = similar(xs)
    v4 = similar(xs)
    safe = similar(xs)
    upstream = similar(xs)

    v1!() = begin
        @inbounds for i in eachindex(xs)
            v1[i] = paper_qw_v1_2pass(xs[i], ys[i], cs[i])
        end
        v1
    end
    v4!() = begin
        @inbounds for i in eachindex(xs)
            v4[i] = paper_qw_v4_5pass(xs[i], ys[i], cs[i])
        end
        v4
    end
    safe!() = begin
        @inbounds for i in eachindex(xs)
            safe[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        safe
    end
    upstream!() = begin
        @inbounds for i in eachindex(xs)
            upstream[i] = xs[i] * ys[i] + cs[i]
        end
        upstream
    end

    v1!(); v4!(); safe!(); upstream!()
    t1 = minimum_time(v1!)
    t4 = minimum_time(v4!)
    ts = minimum_time(safe!)
    tu = minimum_time(upstream!)

    println("QW / Float64x4 scalar timing — $n ordinary cases")
    @printf("  paper v1 former 146-flop / 2-pass = %.3f ms\n", 1e3*t1)
    @printf("  paper v4 proved 176-flop / 5-pass = %.3f ms\n", 1e3*t4)
    @printf("  project safe + renormalize          = %.3f ms\n", 1e3*ts)
    @printf("  upstream MultiFloats x*y+c          = %.3f ms\n", 1e3*tu)
    @printf("  v1/v4 slowdown = %.3fx; safe/v4 = %.3fx; upstream/v4 = %.3fx\n", t4/t1, ts/t4, tu/t4)

    return xs, ys, cs
end

function audit_qw_simd(::Val{W}; scalar_lanes=20_000) where {W}
    Random.seed!(0x2607_5000 + W)
    sx = [signed_rand(T4P) for _ in 1:scalar_lanes]
    sy = [signed_rand(T4P) for _ in 1:scalar_lanes]
    sc = [signed_rand(T4P) for _ in 1:scalar_lanes]
    xs = pack_vectors(Val(W), sx)
    ys = pack_vectors(Val(W), sy)
    cs = pack_vectors(Val(W), sc)
    v1 = similar(xs)
    v4 = similar(xs)
    safe = similar(xs)
    upstream = similar(xs)

    v1!() = begin
        @inbounds for i in eachindex(xs)
            v1[i] = paper_qw_v1_2pass(xs[i], ys[i], cs[i])
        end
        v1
    end
    v4!() = begin
        @inbounds for i in eachindex(xs)
            v4[i] = paper_qw_v4_5pass(xs[i], ys[i], cs[i])
        end
        v4
    end
    safe!() = begin
        @inbounds for i in eachindex(xs)
            safe[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        safe
    end
    upstream!() = begin
        @inbounds for i in eachindex(xs)
            upstream[i] = xs[i] * ys[i] + cs[i]
        end
        upstream
    end

    t1 = minimum_time(v1!)
    t4 = minimum_time(v4!)
    ts = minimum_time(safe!)
    tu = minimum_time(upstream!)
    @printf("  Vec%d: v1=%.3f ms, v4=%.3f ms, safe=%.3f ms, upstream=%.3f ms; upstream/v4=%.3fx\n",
            W, 1e3*t1, 1e3*t4, 1e3*ts, 1e3*tu, tu/t4)
end

function main()
    println("arXiv:2607.11391 v4 reproduction/audit; informational")
    println("CPU: ", Sys.CPU_NAME)
    println("MultiFloats: ", Base.pkgversion(MultiFloats))
    println("Paper distinction: v1 QW=146 flops/2 passes; current v4 QW=176 flops/5 passes.")
    audit_dw()

    xs, ys, cs = audit_qw_scalar_speed()
    audit_qw_corpus("QW ordinary audit", xs, ys, cs; exact_sample=1_000)

    Random.seed!(0x2607_6004)
    nw = 5_000
    xw = [wide_rand(T4P) for _ in 1:nw]
    yw = [wide_rand(T4P) for _ in 1:nw]
    cw = [wide_rand(T4P) for _ in 1:nw]
    audit_qw_corpus("QW wide-exponent audit", xw, yw, cw; exact_sample=500)

    nc = 5_000
    xc = [wide_rand(T4P; emin=-180, emax=180) for _ in 1:nc]
    yc = [wide_rand(T4P; emin=-180, emax=180) for _ in 1:nc]
    cc = Vector{T4P}(undef, nc)
    setprecision(BigFloat, 1024) do
        @inbounds for i in 1:nc
            bits = (20, 50, 80, 110)[mod1(i, 4)]
            delta = BigFloat(2)^(-bits)
            xy = BigFloat(xc[i]) * BigFloat(yc[i])
            cc[i] = T4P(-xy * (one(BigFloat) - delta))
        end
    end
    audit_qw_corpus("QW destructive-cancellation audit", xc, yc, cc; exact_sample=1_000)

    println("QW SIMD timing — 20,000 scalar lanes")
    for W in (2, 4, 8)
        audit_qw_simd(Val(W))
    end
end

main()
