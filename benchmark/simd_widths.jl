using MultiFloatArithmetic
using MultiFloats
using Random

const CASES = (
    (2, Float64x2),
    (3, Float64x3),
    (4, Float64x4),
)

function minimum_time(f; samples=7)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function make_vecs(::Type{T}, ::Val{N}, ::Val{W}, scalar_lanes) where {T,N,W}
    nvec = scalar_lanes ÷ W
    V = MultiFloatVec{W,Float64,N}
    scalars = rand(T, W * nvec)
    vecs = Vector{V}(undef, nvec)
    @inbounds for i in 1:nvec
        base = W * (i - 1)
        vecs[i] = V(ntuple(lane -> scalars[base + lane], W))
    end
    return vecs
end

function bench_width(::Type{T}, ::Val{N}, ::Val{W}; scalar_lanes=48_000) where {T,N,W}
    Random.seed!(0x51d0_2026 + W + 16N)
    xs = make_vecs(T, Val(N), Val(W), scalar_lanes)
    ys = make_vecs(T, Val(N), Val(W), scalar_lanes)
    cs = make_vecs(T, Val(N), Val(W), scalar_lanes)
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

    fused!(); separate!()
    tf = minimum_time(fused!)
    ts = minimum_time(separate!)
    throughput = scalar_lanes / tf / 1e6
    println("$(T) Vec$(W): fused=$(round(tf*1e3; digits=3)) ms, ",
            "mul+add=$(round(ts*1e3; digits=3)) ms, ",
            "speedup=$(round(ts/tf; digits=3))x, ",
            "fused_throughput=$(round(throughput; digits=1)) M scalar-FMA/s")
end

println("SIMD-width smoke benchmark; informational, not a hard performance gate")
println("CPU: ", Sys.CPU_NAME)
for (N, T) in CASES
    for W in (2, 4, 8)
        bench_width(T, Val(N), Val(W))
    end
end
