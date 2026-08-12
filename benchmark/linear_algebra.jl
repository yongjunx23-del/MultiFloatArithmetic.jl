using LinearAlgebra
using MultiFloatArithmetic
using MultiFloats
using Random
using Printf

const MFLA = MultiFloatArithmetic.MFLinearAlgebra

function minimum_time(f; samples=5)
    f()
    best = Inf
    for _ in 1:samples
        GC.gc()
        best = min(best, @elapsed f())
    end
    return best
end

function bench_gemm(::Type{M}, n) where {M<:MultiFloat}
    Random.seed!(0x6e6d + n + sizeof(M))
    A = rand(M, n, n)
    B = rand(M, n, n)
    Cpublic = zeros(M, n, n)
    Cstream = similar(Cpublic)
    Cvec8 = similar(Cpublic)
    Cvec8x2 = similar(Cpublic)
    Cgeneric = similar(Cpublic)

    public!() = MFLA.gemm!(Cpublic, A, B)
    streaming!() = MFLA._gemm_streaming_kernel!(Cstream, A, B)
    vec8!() = MFLA._gemm_vec!(Cvec8, A, B, Val(8))
    vec8x2!() = MFLA._gemm_vec2col!(Cvec8x2, A, B, Val(8))
    generic!() = mul!(Cgeneric, A, B)

    public!(); streaming!(); vec8!(); vec8x2!(); generic!()
    @assert all(MultiFloats.isnormalized, Cpublic)
    @assert all(MultiFloats.isnormalized, Cstream)
    @assert all(MultiFloats.isnormalized, Cvec8)
    @assert all(MultiFloats.isnormalized, Cvec8x2)
    @assert all(i -> Cpublic[i] === Cstream[i], eachindex(Cstream))
    @assert all(i -> Cvec8[i] === Cstream[i], eachindex(Cstream))
    @assert all(i -> Cvec8x2[i] === Cstream[i], eachindex(Cstream))

    tp = minimum_time(public!)
    ts = minimum_time(streaming!)
    t8 = minimum_time(vec8!)
    t82 = minimum_time(vec8x2!)
    tg = minimum_time(generic!)
    work = 2n^3
    @printf("  GEMM n=%d: public=%8.3f ms (%7.1f MFLOP/s), streaming=%8.3f ms (public speedup=%5.2fx), Vec8=%8.3f ms, Vec8x2=%8.3f ms, generic=%8.3f ms, public/generic=%5.2fx\n",
            n, 1e3tp, work/tp/1e6, 1e3ts, ts/tp, 1e3t8, 1e3t82,
            1e3tg, tg/tp)
end

# Exercise row and column remainders independently from the square timing cases.
function check_remainders(::Type{M}) where {M<:MultiFloat}
    Random.seed!(0x8e_2026 + sizeof(M))
    A = rand(M, 19, 11)
    B = rand(M, 11, 7)
    Cstream = zeros(M, 19, 7)
    Cpublic = similar(Cstream)
    MFLA._gemm_streaming_kernel!(Cstream, A, B)
    MFLA.gemm!(Cpublic, A, B)
    @assert all(i -> Cpublic[i] === Cstream[i], eachindex(Cstream))
    @assert all(MultiFloats.isnormalized, Cpublic)
end

function bench_gemv(::Type{M}, m, n) where {M<:MultiFloat}
    Random.seed!(0x6e76 + m + n + sizeof(M))
    A = rand(M, m, n)
    x = rand(M, n)
    yfast = zeros(M, m)
    ygeneric = similar(yfast)

    fast!() = MFLA.gemv!(yfast, A, x)
    generic!() = mul!(ygeneric, A, x)

    fast!(); generic!()
    tf = minimum_time(fast!)
    tg = minimum_time(generic!)
    work = 2m*n
    @printf("  GEMV %dx%d: direct-FMA=%8.3f ms (%7.1f MFLOP/s), generic=%8.3f ms, speedup=%5.2fx\n",
            m, n, 1e3tf, work/tf/1e6, 1e3tg, tg/tf)
end

function bench_dot(::Type{M}, n) where {M<:MultiFloat}
    Random.seed!(0xd07 + n + sizeof(M))
    x = rand(M, n)
    y = rand(M, n)
    sink = Ref(zero(M))

    fast!() = (sink[] = MFLA.mfdot(x, y))
    generic!() = (sink[] = sum(x .* y))

    tf = minimum_time(fast!)
    tg = minimum_time(generic!)
    @printf("  DOT n=%d: direct-FMA=%8.3f ms, generic mul+sum=%8.3f ms, speedup=%5.2fx\n",
            n, 1e3tf, 1e3tg, tg/tf)
end

println("MultiFloat linear algebra direct-FMA benchmark; informational")
println("CPU: ", Sys.CPU_NAME)
println("MultiFloats: ", Base.pkgversion(MultiFloats))
for M in (Float64x2, Float64x4)
    println(M)
    check_remainders(M)
    bench_dot(M, 20_000)
    bench_gemv(M, 256, 128)
    for n in (16, 32, 48, 64)
        bench_gemm(M, n)
    end
end