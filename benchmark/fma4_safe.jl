using MultiFloatArithmetic
using MultiFloats
using Random

import MultiFloats: MultiFloat, MultiFloatVec, fast_two_sum, two_prod, two_sum

# Executable copy of the pre-repair x4 end network. It is retained only to
# measure the cost and numerical effect of the cancellation-safe baseline.
@inline function fma4_old_limbs(
    x::NTuple{4,T},
    y::NTuple{4,T},
    c::NTuple{4,T},
) where {T}
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

    w0, w1 = fast_two_sum(b, a1)
    w1, w2 = two_sum(w1, a2)
    w2, w3 = two_sum(w2, a3)
    z0, rho = two_sum(w0, w1)
    z1, sigma = two_sum(rho, w2)
    z2, z3 = fast_two_sum(sigma, w3)
    return (z0, z1, z2, z3)
end

@inline fma4_old(x::MultiFloat{T,4}, y::MultiFloat{T,4}, c::MultiFloat{T,4}) where {T} =
    MultiFloat{T,4}(fma4_old_limbs(x._limbs, y._limbs, c._limbs))

@inline fma4_old(
    x::MultiFloatVec{W,T,4},
    y::MultiFloatVec{W,T,4},
    c::MultiFloatVec{W,T,4},
) where {W,T} = MultiFloatVec{W,T,4}(
    fma4_old_limbs(x._limbs, y._limbs, c._limbs),
)

# Inspect the exact precondition used by the first FastTwoSum in the old x4
# network. This is diagnostic only and is evaluated on concrete Float64 inputs.
@inline function old_first_fast_precondition(x::Float64x4, y::Float64x4, c::Float64x4)
    x0, x1, _, _ = x._limbs
    y0, y1, _, _ = y._limbs
    c0, c1, _, _ = c._limbs

    p00, e00 = two_prod(x0, y0)
    p01, _ = two_prod(x0, y1)
    p10, _ = two_prod(x1, y0)

    b, r = two_sum(p00, c0)
    a1, _ = two_sum(p01, p10)
    a1, _ = two_sum(a1, e00)
    a1, _ = two_sum(a1, c1)
    a1, _ = two_sum(a1, r)
    return abs(b) >= abs(a1)
end

function minimum_time(f; samples=7)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function benchmark_scalar(; n=20_000)
    Random.seed!(0xf4a0_2026)
    T = Float64x4
    xs = rand(T, n)
    ys = rand(T, n)
    cs = rand(T, n)
    old = similar(xs)
    safe = similar(xs)
    upstream = similar(xs)

    old!() = begin
        @inbounds for i in eachindex(xs)
            old[i] = fma4_old(xs[i], ys[i], cs[i])
        end
        old
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

    old!(); safe!(); upstream!()
    @assert all(MultiFloats.isnormalized, safe)

    to = minimum_time(old!)
    ts = minimum_time(safe!)
    tu = minimum_time(upstream!)
    changed = count(i -> old[i] !== safe[i], eachindex(old))

    println("Float64x4 scalar: old=$(round(to*1e3; digits=3)) ms, ",
            "safe=$(round(ts*1e3; digits=3)) ms, ",
            "upstream=$(round(tu*1e3; digits=3)) ms, ",
            "safe/old=$(round(ts/to; digits=3))x, ",
            "upstream/safe=$(round(tu/ts; digits=3))x, changed=$(changed)/$(n)")
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

function benchmark_vector(::Val{W}; scalar_lanes=20_000) where {W}
    Random.seed!(0xf4a0_2026 + W)
    T = Float64x4
    xs = pack_vectors(Val(W), rand(T, scalar_lanes))
    ys = pack_vectors(Val(W), rand(T, scalar_lanes))
    cs = pack_vectors(Val(W), rand(T, scalar_lanes))
    old = similar(xs)
    safe = similar(xs)

    old!() = begin
        @inbounds for i in eachindex(xs)
            old[i] = fma4_old(xs[i], ys[i], cs[i])
        end
        old
    end
    safe!() = begin
        @inbounds for i in eachindex(xs)
            safe[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        safe
    end

    old!(); safe!()
    changed = 0
    for i in eachindex(safe), lane in 1:W
        @assert MultiFloats.isnormalized(safe[i][lane])
        changed += old[i][lane] !== safe[i][lane]
    end

    to = minimum_time(old!)
    ts = minimum_time(safe!)
    println("Float64x4 Vec$(W): old=$(round(to*1e3; digits=3)) ms, ",
            "safe=$(round(ts*1e3; digits=3)) ms, safe/old=$(round(ts/to; digits=3))x, ",
            "changed=$(changed)/$(scalar_lanes)")
end

function cancellation_diagnostic(; n=10_000)
    Random.seed!(0xc4nc_2026)
    T = Float64x4
    violations = 0
    changed = 0
    old_bad = 0
    safe_bad = 0

    setprecision(BigFloat, 1024) do
        for _ in 1:n
            x = rand(T)
            y = rand(T)
            c = -T(big(x) * big(y))
            violations += !old_first_fast_precondition(x, y, c)

            old = fma4_old(x, y, c)
            safe = fma_fast(x, y, c)
            changed += old !== safe
            old_bad += !MultiFloats.isnormalized(old)
            safe_bad += !MultiFloats.isnormalized(safe)
        end
    end

    @assert safe_bad == 0
    println("cancellation diagnostic: first-FastTwoSum violations=$(violations)/$(n), ",
            "old_vs_safe_changed=$(changed)/$(n), old_non_normalized=$(old_bad), ",
            "safe_non_normalized=$(safe_bad)")
end

println("Cancellation-safe Float64x4 baseline benchmark; informational")
println("CPU: ", Sys.CPU_NAME)
benchmark_scalar()
for W in (2, 4, 8)
    benchmark_vector(Val(W))
end
cancellation_diagnostic()
