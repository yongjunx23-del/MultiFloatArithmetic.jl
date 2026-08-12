using MultiFloatArithmetic
using MultiFloats
using Random

const ExperimentalArithmetic = MultiFloatArithmetic.Experimental

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

@inline one_limb_expansion(q::T, ::Val{N}) where {T,N} =
    ntuple(i -> isone(i) ? q : zero(T), Val{N}())

@inline function div_digits_full_limbs(x::NTuple{N,T}, y::NTuple{N,T}) where {N,T}
    residual = x
    digits = ntuple(_ -> zero(T), Val{N}())
    @inbounds for k in 1:N
        qk = MultiFloats.div_r(residual[1], y[1])
        digits = Base.setindex(digits, qk, k)
        qexp = one_limb_expansion(qk, Val{N}())
        product = MultiFloats.mfmul(y, qexp, Val{N}())
        residual = MultiFloats.mfadd(residual, map(-, product), Val{N}())
    end
    return MultiFloats.renormalize(digits)
end

@inline div_digits_full(x::MultiFloat{T,N}, y::MultiFloat{T,N}) where {T,N} =
    MultiFloat{T,N}(div_digits_full_limbs(x._limbs, y._limbs))

@inline div_digits_full(x::MultiFloatVec{W,T,N}, y::MultiFloatVec{W,T,N}) where {W,T,N} =
    MultiFloatVec{W,T,N}(div_digits_full_limbs(x._limbs, y._limbs))

function report_division(label, specialized!, full!, upstream!)
    specialized!(); full!(); upstream!()
    ts = min_time(specialized!)
    tf = min_time(full!)
    tu = min_time(upstream!)
    println("$(label): specialized=$(round(ts*1e3; digits=3)) ms, ",
            "full=$(round(tf*1e3; digits=3)) ms, ",
            "upstream=$(round(tu*1e3; digits=3)) ms, ",
            "full/specialized=$(round(tf/ts; digits=3))x, ",
            "upstream/specialized=$(round(tu/ts; digits=3))x")
end

function bench_scalar(::Type{T}; n=12_000) where {T}
    Random.seed!(0xd1a1_2026)
    xs = rand(T, n) .+ T(0.5)
    ys = rand(T, n) .+ T(0.5)
    out = similar(xs)

    specialized!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = ExperimentalArithmetic.div_digits(xs[i], ys[i])
        end
        out
    end

    full!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = div_digits_full(xs[i], ys[i])
        end
        out
    end

    upstream!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = xs[i] / ys[i]
        end
        out
    end

    report_division("$(T) scalar", specialized!, full!, upstream!)
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

    specialized!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = ExperimentalArithmetic.div_digits(xs[i], ys[i])
        end
        out
    end

    full!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = div_digits_full(xs[i], ys[i])
        end
        out
    end

    upstream!() = begin
        @inbounds for i in eachindex(xs)
            out[i] = xs[i] / ys[i]
        end
        out
    end

    report_division("$(T) Vec4", specialized!, full!, upstream!)
end

println("Rejected quotient-digit division A/B; informational and retained for reproducibility")
println("CPU: ", Sys.CPU_NAME)
for (N, T) in DIV_CASES
    bench_scalar(T)
    bench_vec4(T, Val(N))
end
