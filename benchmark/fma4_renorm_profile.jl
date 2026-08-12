using MultiFloatArithmetic
using MultiFloats
using Random
using Printf

import MultiFloats: two_prod, two_sum, renormalize

const T4P = MultiFloats.Float64x4

# Reproduce the accepted x4 direct-FMA network only up to the four-term
# pre-renormalization tuple. This is a diagnostic copy: production remains in
# src/MultiFloatArithmetic.jl.
@inline function fma4_pre_limbs(x::NTuple{4,T}, y::NTuple{4,T}, c::NTuple{4,T}) where {T}
    x0, x1, x2, x3 = x
    y0, y1, y2, y3 = y
    c0, c1, c2, c3 = c

    p00, e00 = two_prod(x0, y0)
    p01, e01 = two_prod(x0, y1)
    p10, e10 = two_prod(x1, y0)
    p02, e02 = two_prod(x0, y2)
    p11, e11 = two_prod(x1, y1)
    p20, e20 = two_prod(x2, y0)

    p03 = x0 * y3
    p12 = x1 * y2
    p21 = x2 * y1
    p30 = x3 * y0
    diagonal3 = (p03 + p30) + (p12 + p21)

    b, r = two_sum(p00, c0)

    a1, f1 = two_sum(p01, p10)
    a1, f2 = two_sum(a1, e00)
    a1, f3 = two_sum(a1, c1)
    a1, f4 = two_sum(a1, r)

    a2, g1 = two_sum(p02, p20)
    a2, g2 = two_sum(a2, p11)
    e01e10, g4 = two_sum(e01, e10)
    a2, g3 = two_sum(a2, e01e10)
    a2, g5 = two_sum(a2, c2)
    a2, g6 = two_sum(a2, f1)
    a2, g7 = two_sum(a2, f2)
    a2, g8 = two_sum(a2, f3)
    a2, g9 = two_sum(a2, f4)

    t1 = e02 + e20
    t2 = e11 + diagonal3
    t3 = (t1 + t2) + c3

    t1 = g1 + g2
    t2 = g3 + g4
    t1 = t1 + t2
    t2 = g6 + g7
    t4 = g8 + g9
    t2 = t2 + t4
    t1 = t1 + t2
    t1 = t1 + g5
    a3 = t3 + t1

    w0, w1 = two_sum(b, a1)
    w1, w2 = two_sum(w1, a2)
    w2, w3 = two_sum(w2, a3)
    z0, rho = two_sum(w0, w1)
    z1, sigma = two_sum(rho, w2)
    z2, z3 = two_sum(sigma, w3)
    return (z0, z1, z2, z3)
end

# Exact N=4 equivalent of MultiFloats._renorm_pass, kept local so this
# benchmark does not depend on an upstream private symbol.
@inline function renorm_pass4(x::NTuple{4,T}) where {T}
    x1, x2, x3, x4 = x
    x1, x2 = two_sum(x1, x2)
    x3, x4 = two_sum(x3, x4)
    x2, x3 = two_sum(x2, x3)
    return (x1, x2, x3, x4)
end

function renorm_update_count(x::NTuple{4,T}; max_updates=32) where {T}
    current = x
    for updates in 0:max_updates
        next = renorm_pass4(current)
        next === current && return updates
        current = next
    end
    return max_updates + 1
end

# Candidate shape only: unroll the common first three changing passes, then
# retain the authoritative generic fallback. This benchmark establishes whether
# such a specialization has enough coverage and speed potential to justify a
# separately verified production implementation.
@inline function renorm4_unrolled3_fallback(x::NTuple{4,T}) where {T}
    p1 = renorm_pass4(x)
    p1 === x && return p1
    p2 = renorm_pass4(p1)
    p2 === p1 && return p2
    p3 = renorm_pass4(p2)
    p3 === p2 && return p3
    return renormalize(p3)
end

@inline fma4_unrolled3(x::T4P, y::T4P, c::T4P) =
    T4P(renorm4_unrolled3_fallback(fma4_pre_limbs(x._limbs, y._limbs, c._limbs)))

function minimum_time(f; samples=7)
    f()
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function profile_updates(label, xs, ys, cs)
    counts = Dict{Int,Int}()
    max_seen = 0
    for i in eachindex(xs)
        pre = fma4_pre_limbs(xs[i]._limbs, ys[i]._limbs, cs[i]._limbs)
        n = renorm_update_count(pre)
        counts[n] = get(counts, n, 0) + 1
        max_seen = max(max_seen, n)
    end
    total = length(xs)
    covered3 = sum(v for (k, v) in counts if k <= 3)
    println(label)
    println("  cases=", total, ", max changing passes=", max_seen,
            ", covered by <=3 changing passes=", covered3, "/", total)
    println("  histogram changing_passes => cases: ", sort(collect(counts)))
    return counts
end

function cancellation_accuracy(xs, ys, cs)
    direct_better = 0
    upstream_better = 0
    equal_error = 0
    bitwise_same = 0
    direct_exact = 0
    upstream_exact = 0
    for i in eachindex(xs)
        direct = fma_fast(xs[i], ys[i], cs[i])
        upstream = xs[i] * ys[i] + cs[i]
        ref = Rational{BigInt}(xs[i]) * Rational{BigInt}(ys[i]) + Rational{BigInt}(cs[i])
        ed = abs(Rational{BigInt}(direct) - ref)
        eu = abs(Rational{BigInt}(upstream) - ref)
        direct_better += ed < eu
        upstream_better += eu < ed
        equal_error += ed == eu
        bitwise_same += direct === upstream
        direct_exact += iszero(ed)
        upstream_exact += iszero(eu)
    end
    println("destructive-cancellation direct-FMA vs upstream x*y+c")
    println("  direct_better=", direct_better, "/", length(xs),
            ", upstream_better=", upstream_better,
            ", equal_error=", equal_error)
    println("  bitwise_same=", bitwise_same, "/", length(xs),
            ", direct_exact=", direct_exact,
            ", upstream_exact=", upstream_exact)
end

function main()
    Random.seed!(0xf4a4_2026)
    n = 20_000
    xs = rand(T4P, n)
    ys = rand(T4P, n)
    cs = rand(T4P, n)

    profile_updates("ordinary x4 final-renormalization profile", xs, ys, cs)

    nc = 5_000
    xcancel = rand(T4P, nc)
    ycancel = rand(T4P, nc)
    ccancel = Vector{T4P}(undef, nc)
    setprecision(BigFloat, 1024) do
        @inbounds for i in eachindex(ccancel)
            ccancel[i] = T4P(-big(xcancel[i]) * big(ycancel[i]))
        end
    end
    profile_updates("destructive-cancellation x4 profile", xcancel, ycancel, ccancel)
    cancellation_accuracy(xcancel, ycancel, ccancel)

    current = similar(xs)
    specialized = similar(xs)
    upstream = similar(xs)
    current!() = begin
        @inbounds for i in eachindex(xs)
            current[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        current
    end
    specialized!() = begin
        @inbounds for i in eachindex(xs)
            specialized[i] = fma4_unrolled3(xs[i], ys[i], cs[i])
        end
        specialized
    end
    upstream!() = begin
        @inbounds for i in eachindex(xs)
            upstream[i] = xs[i] * ys[i] + cs[i]
        end
        upstream
    end

    current!(); specialized!(); upstream!()
    @assert specialized == current

    # Cancellation equality is the non-negotiable safety check for this
    # diagnostic specialization.
    @inbounds for i in eachindex(xcancel)
        @assert fma4_unrolled3(xcancel[i], ycancel[i], ccancel[i]) ===
                fma_fast(xcancel[i], ycancel[i], ccancel[i])
    end

    tc = minimum_time(current!)
    ts = minimum_time(specialized!)
    tu = minimum_time(upstream!)

    println("Float64x4 scalar finalizer A/B (", n, " ordinary cases)")
    @printf("  current safe fma_fast = %.3f ms\n", 1e3 * tc)
    @printf("  unrolled3 + fallback = %.3f ms\n", 1e3 * ts)
    @printf("  upstream x*y+c       = %.3f ms\n", 1e3 * tu)
    @printf("  current/specialized   = %.3fx speedup potential\n", tc / ts)
    @printf("  specialized/upstream  = %.3fx time ratio\n", ts / tu)
end

println("Float64x4 renormalization specialization diagnostic; informational")
println("CPU: ", Sys.CPU_NAME)
main()
