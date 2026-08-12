using MultiFloatArithmetic
using MultiFloats
using Random

import MultiFloats: MultiFloat, MultiFloatVec, fast_two_sum, two_prod, two_sum

# Executable copy of the pre-repair x3 end network. This is benchmark-only and
# exists to measure the exact fixed cost of the proof-driven normalization fix.
@inline function fma3_unrepaired_limbs(
    x::NTuple{3,T},
    y::NTuple{3,T},
    c::NTuple{3,T},
) where {T}
    x0, x1, x2 = x
    y0, y1, y2 = y
    c0, c1, c2 = c

    p00, e00 = two_prod(x0, y0)
    p01, e01 = two_prod(x0, y1)
    p10, e10 = two_prod(x1, y0)

    p02 = x0 * y2
    p11 = x1 * y1
    p20 = x2 * y0

    sigma = (p02 + p20) + p11
    tail = ((e01 + e10) + sigma) + c2

    a, q1 = two_sum(p01, p10)
    a, q2 = two_sum(a, e00)
    a, q3 = two_sum(a, c1)
    tail = tail + ((q1 + q2) + q3)

    b, r = two_sum(p00, c0)
    m1, m2 = two_sum(r, a)
    m2 = m2 + tail

    w0, w1 = fast_two_sum(b, m1)
    w1, w2 = two_sum(w1, m2)
    z0, rho = two_sum(w0, w1)
    z1, z2 = fast_two_sum(rho, w2)
    return (z0, z1, z2)
end

@inline fma3_unrepaired(
    x::MultiFloat{T,3},
    y::MultiFloat{T,3},
    c::MultiFloat{T,3},
) where {T} = MultiFloat{T,3}(fma3_unrepaired_limbs(x._limbs, y._limbs, c._limbs))

@inline fma3_unrepaired(
    x::MultiFloatVec{W,T,3},
    y::MultiFloatVec{W,T,3},
    c::MultiFloatVec{W,T,3},
) where {W,T} = MultiFloatVec{W,T,3}(
    fma3_unrepaired_limbs(x._limbs, y._limbs, c._limbs),
)

function minimum_time(f; samples=7)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function validate_scalar(old, repaired)
    setprecision(BigFloat, 512) do
        for i in eachindex(old)
            @assert big(old[i]) == big(repaired[i])
            @assert MultiFloats.isnormalized(repaired[i])
        end
    end
end

function report(label, old!, repaired!, old, repaired)
    old!()
    repaired!()
    validate_scalar(old, repaired)

    to = minimum_time(old!)
    tr = minimum_time(repaired!)
    old_bad = count(x -> !MultiFloats.isnormalized(x), old)
    repaired_bad = count(x -> !MultiFloats.isnormalized(x), repaired)

    println("$(label): old=$(round(to * 1e3; digits=3)) ms, ",
            "repaired=$(round(tr * 1e3; digits=3)) ms, ",
            "repair/old=$(round(tr / to; digits=3))x, ",
            "observed_non_normalized=$(old_bad)/$(repaired_bad) old/repaired")
end

function benchmark_scalar(; n=20_000)
    Random.seed!(0xf3a0_2026)
    T = Float64x3
    xs = rand(T, n)
    ys = rand(T, n)
    cs = rand(T, n)
    old = similar(xs)
    repaired = similar(xs)

    old!() = begin
        @inbounds for i in eachindex(xs)
            old[i] = fma3_unrepaired(xs[i], ys[i], cs[i])
        end
        old
    end

    repaired!() = begin
        @inbounds for i in eachindex(xs)
            repaired[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        repaired
    end

    report("Float64x3 scalar", old!, repaired!, old, repaired)
end

function pack_vectors(::Val{W}, scalars) where {W}
    V = MultiFloatVec{W,Float64,3}
    nvec = length(scalars) ÷ W
    out = Vector{V}(undef, nvec)
    @inbounds for i in 1:nvec
        base = W * (i - 1)
        out[i] = V(ntuple(lane -> scalars[base + lane], Val(W)))
    end
    return out
end

function validate_vectors(old, repaired, ::Val{W}) where {W}
    setprecision(BigFloat, 512) do
        for i in eachindex(old), lane in 1:W
            @assert big(old[i][lane]) == big(repaired[i][lane])
            @assert MultiFloats.isnormalized(repaired[i][lane])
        end
    end
end

function benchmark_vector(::Val{W}; scalar_lanes=20_000) where {W}
    Random.seed!(0xf3a0_2026 + W)
    T = Float64x3
    xs = pack_vectors(Val(W), rand(T, scalar_lanes))
    ys = pack_vectors(Val(W), rand(T, scalar_lanes))
    cs = pack_vectors(Val(W), rand(T, scalar_lanes))
    old = similar(xs)
    repaired = similar(xs)

    old!() = begin
        @inbounds for i in eachindex(xs)
            old[i] = fma3_unrepaired(xs[i], ys[i], cs[i])
        end
        old
    end

    repaired!() = begin
        @inbounds for i in eachindex(xs)
            repaired[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        repaired
    end

    old!()
    repaired!()
    validate_vectors(old, repaired, Val(W))

    to = minimum_time(old!)
    tr = minimum_time(repaired!)
    old_bad = 0
    repaired_bad = 0
    for i in eachindex(old), lane in 1:W
        old_bad += !MultiFloats.isnormalized(old[i][lane])
        repaired_bad += !MultiFloats.isnormalized(repaired[i][lane])
    end

    println("Float64x3 Vec$(W): old=$(round(to * 1e3; digits=3)) ms, ",
            "repaired=$(round(tr * 1e3; digits=3)) ms, ",
            "repair/old=$(round(tr / to; digits=3))x, ",
            "observed_non_normalized=$(old_bad)/$(repaired_bad) old/repaired")
end

println("Proof-driven Float64x3 final-normalization A/B; informational")
println("CPU: ", Sys.CPU_NAME)
benchmark_scalar()
for W in (2, 4, 8)
    benchmark_vector(Val(W))
end
