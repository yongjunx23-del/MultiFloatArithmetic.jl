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
    Cfast = zeros(M, n, n)
    Cvec4 = similar(Cfast)
    Cvec8 = similar(Cfast)
    Cvec8x2 = similar(Cfast)
    Cgeneric = similar(Cfast)

    fast!() = MFLA.gemm!(Cfast, A, B)
    vec4!() = MFLA._gemm_vec!(Cvec4, A, B, Val(4))
    vec8!() = MFLA._gemm_vec!(Cvec8, A, B, Val(8))
    vec8x2!() = MFLA._gemm_vec2col!(Cvec8x2, A, B, Val(8))
    generic!() = mul!(Cgeneric, A, B)

    fast!(); vec4!(); vec8!(); vec8x2!(); generic!()
    @assert all(MultiFloats.isnormalized, Cfast)
    @assert all(MultiFloats.isnormalized, Cvec4)
    @assert all(MultiFloats.isnormalized, Cvec8)
    @assert all(MultiFloats.isnormalized, Cvec8x2)
    @assert all(i -> Cvec4[i] === Cfast[i], eachindex(Cfast))
    @assert all(i -> Cvec8[i] === Cfast[i], eachindex(Cfast))
    @assert all(i -> Cvec8x2[i] === Cfast[i], eachindex(Cfast))

    tf = minimum_time(fast!)
    t4 = minimum_time(vec4!)
    t8 = minimum_time(vec8!)
    t82 = minimum_time(vec8x2!)
    tg = minimum_time(generic!)
    work = 2n^3
    @printf("  GEMM n=%d: streaming=%8.3f ms (%7.1f MFLOP/s), Vec4=%8.3f ms (%5.2fx), Vec8=%8.3f ms (%5.2fx), Vec8x2=%8.3f ms (%5.2fx vs streaming, %5.2fx vs Vec8), generic=%8.3f ms, streaming/generic=%5.2fx\n",
            n, 1e3tf, work/tf/1e6, 1e3t4, tf/t4, 1e3t8, tf/t8,
            1e3t82, tf/t82, t8/t82, 1e3tg, tg/tf)
end

# Exercise row and column remainders independently from the square timing cases.
function check_remainders(::Type{M}) where {M<:MultiFloat}
    Random.seed!(0x8e_2026 + sizeof(M))
    A = rand(M, 19, 11)
    B = rand(M, 11, 7)
    Cstream = zeros(M, 19, 7)
    Cvec = similar(Cstream)
    MFLA.gemm!(Cstream, A, B)
    MFLA._gemm_vec2col!(Cvec, A, B, Val(8))
    @assert all(i -> Cvec[i] === Cstream[i], eachindex(Cstream))
    @assert all(MultiFloats.isnormalized, Cvec)
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