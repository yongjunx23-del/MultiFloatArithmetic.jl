using MultiFloatArithmetic
using MultiFloats
using Random
using Printf

include(joinpath(@__DIR__, "..", "audit", "paper2607_v1.jl"))
using .Paper2607V1Audit

const T2A = Float64x2
const T4A = Float64x4

function minimum_time(f; samples=7)
    f()
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

exact_value(x) = Rational{BigInt}(x)

function signed_rand(::Type{T}) where {T}
    x = rand(T)
    return MultiFloats.renormalize(rand(Bool) ? x : -x)
end

function wide_rand(::Type{T}; emin=-350, emax=350) where {T}
    MultiFloats.renormalize(ldexp(signed_rand(T), rand(emin:emax)))
end

function bound_ratio(z, x, y, c, K)
    truth = exact_value(x) * exact_value(y) + exact_value(c)
    scale = abs(exact_value(x) * exact_value(y)) + abs(exact_value(c))
    err = abs(exact_value(z) - truth)
    iszero(scale) && return BigFloat(0)
    return setprecision(BigFloat, 1024) do
        u = BigFloat(2)^(-53)
        BigFloat(err) / (u^K * BigFloat(scale))
    end
end

function audit_dw(; n=20_000)
    Random.seed!(0x2607_d002)
    xs = [signed_rand(T2A) for _ in 1:n]
    ys = [signed_rand(T2A) for _ in 1:n]
    cs = [signed_rand(T2A) for _ in 1:n]
    paper = similar(xs)
    current = similar(xs)
    upstream = similar(xs)

    paper!() = begin
        @inbounds for i in eachindex(xs)
            paper[i] = paper_dw_fma(xs[i], ys[i], cs[i])
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

    println("DW / Float64x2 — official Algorithm 1")
    println("  paper vs project fma_fast bitwise mismatches = $mismatches/$n")
    @printf("  paper 17-flop FMA = %.3f ms / %d\n", 1e3*tp, n)
    @printf("  current fma_fast   = %.3f ms / %d\n", 1e3*tc, n)
    @printf("  upstream x*y+c     = %.3f ms / %d\n", 1e3*tu, n)
    @printf("  upstream/paper time ratio = %.3fx\n", tu/tp)
end

function audit_qw_corpus(label, xs, ys, cs; exact_sample=1_000)
    n = length(xs)
    nonnormalized = 0
    safe_nonnormalized = 0
    changed = 0
    fast1_fail = 0
    fast2_fail = 0
    exact_sum_mismatch = 0
    symmetry_mismatch = 0
    max_C = BigFloat(0)
    direct_better = 0
    upstream_better = 0
    equal_error = 0

    setprecision(BigFloat, 1024) do
        @inbounds for i in eachindex(xs)
            x, y, c = xs[i], ys[i], cs[i]
            paper = paper_qw_fma(x, y, c)
            safe = fma_fast(x, y, c)
            upstream = x * y + c
            f1, f2 = paper_qw_fast_preconditions(x, y, c)

            nonnormalized += !MultiFloats.isnormalized(paper)
            safe_nonnormalized += !MultiFloats.isnormalized(safe)
            changed += paper !== safe
            fast1_fail += !f1
            fast2_fail += !f2
            symmetry_mismatch += paper !== paper_qw_fma(y, x, c)
            max_C = max(max_C, bound_ratio(paper, x, y, c, 4))

            if i <= exact_sample
                rp = exact_value(paper)
                rs = exact_value(safe)
                exact_sum_mismatch += rp != rs
                truth = exact_value(x) * exact_value(y) + exact_value(c)
                ep = abs(rp - truth)
                eu = abs(exact_value(upstream) - truth)
                if ep < eu
                    direct_better += 1
                elseif eu < ep
                    upstream_better += 1
                else
                    equal_error += 1
                end
            end
        end
    end

    println(label)
    println("  cases=$n")
    println("  paper non-normalized=$nonnormalized; current-safe non-normalized=$safe_nonnormalized; paper_vs_safe_changed=$changed")
    println("  paper FastTwoSum exponent-precondition failures: first=$fast1_fail, tail=$fast2_fail")
    println("  x/y symmetry mismatches=$symmetry_mismatch")
    println("  exact component-sum paper-vs-safe mismatches (first $exact_sample)=$exact_sum_mismatch")
    println("  max observed |err|/(u^4*(|xy|+|c|)) = $max_C (paper certified bound: 812)")
    println("  exact-error vs upstream (first $exact_sample): paper_better=$direct_better, upstream_better=$upstream_better, equal=$equal_error")
end

function audit_qw_scalar_speed(; n=20_000)
    Random.seed!(0x2607_4004)
    xs = [signed_rand(T4A) for _ in 1:n]
    ys = [signed_rand(T4A) for _ in 1:n]
    cs = [signed_rand(T4A) for _ in 1:n]
    paper = similar(xs)
    safe = similar(xs)
    upstream = similar(xs)

    paper!() = begin
        @inbounds for i in eachindex(xs)
            paper[i] = paper_qw_fma(xs[i], ys[i], cs[i])
        end
        paper
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

    tp = minimum_time(paper!)
    ts = minimum_time(safe!)
    tu = minimum_time(upstream!)
    @printf("QW scalar / %d ordinary cases: paper=%.3f ms, current-safe=%.3f ms, upstream=%.3f ms; upstream/paper=%.3fx, safe/paper=%.3fx\n",
            n, 1e3*tp, 1e3*ts, 1e3*tu, tu/tp, ts/tp)
    return xs, ys, cs
end

function pack_vectors(::Val{W}, scalars) where {W}
    V = MultiFloatVec{W,Float64,4}
    nvec = length(scalars) ÷ W
    out = Vector{V}(undef, nvec)
    @inbounds for i in 1:nvec
        base = W * (i - 1)
        out[i] = V(ntuple(lane -> scalars[base + lane], Val(W)))
    end
    out
end

function audit_qw_simd(::Val{W}; scalar_lanes=20_000) where {W}
    Random.seed!(0x2607_5000 + W)
    sx = [signed_rand(T4A) for _ in 1:scalar_lanes]
    sy = [signed_rand(T4A) for _ in 1:scalar_lanes]
    sc = [signed_rand(T4A) for _ in 1:scalar_lanes]
    xs = pack_vectors(Val(W), sx)
    ys = pack_vectors(Val(W), sy)
    cs = pack_vectors(Val(W), sc)
    paper = similar(xs)
    safe = similar(xs)
    upstream = similar(xs)

    paper!() = begin
        @inbounds for i in eachindex(xs)
            paper[i] = paper_qw_fma(xs[i], ys[i], cs[i])
        end
        paper
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

    tp = minimum_time(paper!)
    ts = minimum_time(safe!)
    tu = minimum_time(upstream!)
    @printf("  Vec%d / %d lanes: paper=%.3f ms, safe=%.3f ms, upstream=%.3f ms; upstream/paper=%.3fx\n",
            W, scalar_lanes, 1e3*tp, 1e3*ts, 1e3*tu, tu/tp)
end

function main()
    println("Official arXiv:2607.11391v1 DW/QW reproduction audit")
    println("CPU: ", Sys.CPU_NAME)
    println("MultiFloats version: ", Base.pkgversion(MultiFloats))
    println("Paper contract: DW/QW 17/146 flops; certified C2=34, C4=812; QW leading non-overlap is not guaranteed under extreme cancellation.")

    audit_dw()
    xs, ys, cs = audit_qw_scalar_speed()
    audit_qw_corpus("QW ordinary", xs, ys, cs; exact_sample=1_000)

    Random.seed!(0x2607_6004)
    nw = 5_000
    xw = [wide_rand(T4A) for _ in 1:nw]
    yw = [wide_rand(T4A) for _ in 1:nw]
    cw = [wide_rand(T4A) for _ in 1:nw]
    audit_qw_corpus("QW wide exponent", xw, yw, cw; exact_sample=500)

    nc = 5_000
    xc = [wide_rand(T4A; emin=-180, emax=180) for _ in 1:nc]
    yc = [wide_rand(T4A; emin=-180, emax=180) for _ in 1:nc]
    cc = Vector{T4A}(undef, nc)
    setprecision(BigFloat, 1024) do
        @inbounds for i in 1:nc
            bits = (20, 50, 80, 110, 150)[mod1(i, 5)]
            delta = BigFloat(2)^(-bits)
            xy = BigFloat(xc[i]) * BigFloat(yc[i])
            cc[i] = T4A(-xy * (one(BigFloat) - delta))
        end
    end
    audit_qw_corpus("QW destructive cancellation", xc, yc, cc; exact_sample=1_000)

    println("QW SIMD timing")
    for W in (2, 4, 8)
        audit_qw_simd(Val(W))
    end
end

main()
