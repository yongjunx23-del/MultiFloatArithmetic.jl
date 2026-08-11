using MultiFloatArithmetic
using MultiFloats
using Random

const DIV_CASES = (
    (2, Float64x2),
    (3, Float64x3),
    (4, Float64x4),
)

function min_time(f; samples=7)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function bench_scalar(::Type{T}; n=12_000) where {T}
    Random.seed!(0xd1a1_2026)
    xs = rand(T, n) .+ T(0.5)
    ys = rand(T, n) .+ T(0.5)
    out = similar(xs)

    candidate!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = div_digits(xs[i], ys[i])
        end
        out
    end

    upstream!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = xs[i] / ys[i]
        end
        out
    end

    candidate!(); upstream!()
    tc = min_time(candidate!)
    tu = min_time(upstream!)
    println("$(T) scalar: digits=$(round(tc*1e3; digits=3)) ms, ",
            "upstream=$(round(tu*1e3; digits=3)) ms, ",
            "upstream/digits=$(round(tu/tc; digits=3))x")
end

function pack_vec4(::Type{T}, ::Val{N}, scalars) where {T,N}
    V = MultiFloatVec{4,Float64,N}
    nvec = length(scalars) ÷ 4
    vecs = Vector{V}(undef, nvec)
    @inbounds for i in 1:nvec
        b = 4 * (i - 1)
        vecs[i] = V((scalars[b+1], scalars[b+2], scalars[b+3], scalars[b+4]))
    end
    return vecs
end

function bench_vec4(::Type{T}, ::Val{N}; scalar_lanes=12_000) where {T,N}
    Random.seed!(0xd1a1_2026)
    xs = pack_vec4(T, Val(N), rand(T, scalar_lanes) .+ T(0.5))
    ys = pack_vec4(T, Val(N), rand(T, scalar_lanes) .+ T(0.5))
    out = similar(xs)

    candidate!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = div_digits(xs[i], ys[i])
        end
        out
    end

    upstream!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = xs[i] / ys[i]
        end
        out
    end

    candidate!(); upstream!()
    tc = min_time(candidate!)
    tu = min_time(upstream!)
    println("$(T) Vec4: digits=$(round(tc*1e3; digits=3)) ms, ",
            "upstream=$(round(tu*1e3; digits=3)) ms, ",
            "upstream/digits=$(round(tu/tc; digits=3))x")
end

println("Quotient-digit division research A/B; informational")
println("CPU: ", Sys.CPU_NAME)
for (N, T) in DIV_CASES
    bench_scalar(T)
    bench_vec4(T, Val(N))
end
