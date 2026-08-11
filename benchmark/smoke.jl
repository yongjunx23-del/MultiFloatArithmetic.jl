using MultiFloatArithmetic
using MultiFloats
using Random

const TYPES = (
    MultiFloats.Float64x2,
    MultiFloats.Float64x3,
    MultiFloats.Float64x4,
)

function minimum_time(f; samples=5)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function benchmark_type(::Type{T}; n=20_000) where {T}
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

    fused!()
    separate!()

    tf = minimum_time(fused!)
    ts = minimum_time(separate!)

    println("$(T): fused=$(round(tf * 1e3; digits=3)) ms, ",
            "mul+add=$(round(ts * 1e3; digits=3)) ms, ",
            "mul+add/fused=$(round(ts / tf; digits=3))x")
end

println("Hosted-runner smoke benchmark; informational, not a hard performance gate")
for T in TYPES
    benchmark_type(T)
end
