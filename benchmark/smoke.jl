using MultiFloatArithmetic
using MultiFloats
using Random

const CASES = (
    (2, MultiFloats.Float64x2),
    (3, MultiFloats.Float64x3),
    (4, MultiFloats.Float64x4),
)

function minimum_time(f; samples=7)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function report(label, fused!, separate!)
    fused!()
    separate!()
    tf = minimum_time(fused!)
    ts = minimum_time(separate!)
    println("$(label): fused=$(round(tf * 1e3; digits=3)) ms, ",
            "mul+add=$(round(ts * 1e3; digits=3)) ms, ",
            "mul+add/fused=$(round(ts / tf; digits=3))x")
end

function benchmark_scalar(::Type{T}; n=20_000) where {T}
    Random.seed!(0x5d9a_2026)
    xs = rand(T, n)
    ys = rand(T, n)
    cs = rand(T, n)
    out = similar(xs)

    fused!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        out
    end

    separate!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = xs[i] * ys[i] + cs[i]
        end
        out
    end

    report("$(T) scalar", fused!, separate!)
end

function benchmark_vec4(::Type{T}, ::Val{N}; scalar_count=20_000) where {T,N}
    Random.seed!(0x5d9a_2026)
    nvec = scalar_count ÷ 4
    V = MultiFloatVec{4,Float64,N}

    function make_vecs()
        scalars = rand(T, 4 * nvec)
        vecs = Vector{V}(undef, nvec)
        @inbounds for i in 1:nvec
            base = 4 * (i - 1)
            vecs[i] = V((scalars[base + 1], scalars[base + 2],
                         scalars[base + 3], scalars[base + 4]))
        end
        return vecs
    end

    xs = make_vecs()
    ys = make_vecs()
    cs = make_vecs()
    out = similar(xs)

    fused!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = fma_fast(xs[i], ys[i], cs[i])
        end
        out
    end

    separate!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = xs[i] * ys[i] + cs[i]
        end
        out
    end

    report("$(T) Vec4 ($(4 * nvec) scalar lanes)", fused!, separate!)
end

println("Hosted-runner smoke benchmark; informational, not a hard performance gate")
for (N, T) in CASES
    benchmark_scalar(T)
    benchmark_vec4(T, Val(N))
end
