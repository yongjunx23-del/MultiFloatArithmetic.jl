using MultiFloatArithmetic
using MultiFloats
using Random

function _min_time_2019(f; samples=9)
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function bench_prod33(; n=40_000)
    Random.seed!(0x3a19_2019)
    xs = rand(Float64x3, n) .+ Float64x3(0.5)
    ys = rand(Float64x3, n) .+ Float64x3(0.5)
    out = similar(xs)

    literature!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = tw_prod33_fast(xs[i], ys[i])
        end
        out
    end

    upstream!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = xs[i] * ys[i]
        end
        out
    end

    literature!(); upstream!()
    tl = _min_time_2019(literature!)
    tu = _min_time_2019(upstream!)
    normalized = count(MultiFloats.isnormalized, out)
    println("3Prod_fast(3,3): literature=$(round(tl*1e3; digits=3)) ms, ",
            "MultiFloats=$(round(tu*1e3; digits=3)) ms, ",
            "MultiFloats/literature=$(round(tu/tl; digits=3))x, ",
            "normalized_after_upstream=$(normalized)/$(n)")

    literature!()
    normalized = count(MultiFloats.isnormalized, out)
    println("3Prod_fast(3,3) strong-normalized outputs: $(normalized)/$(n)")
end

println("2019 triple-word scalar multiplication A/B; informational")
println("CPU: ", Sys.CPU_NAME)
bench_prod33()
