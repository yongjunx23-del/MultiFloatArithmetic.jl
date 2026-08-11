using MultiFloatArithmetic
using MultiFloats
using Random

function paper_p_nonoverlap(limbs::NTuple{3,Float64})
    @inbounds for i in 1:2
        hi = limbs[i]
        lo = limbs[i + 1]
        if iszero(hi)
            iszero(lo) || return false
        elseif !(abs(lo) < eps(abs(hi)))
            return false
        end
    end
    return true
end

function wide_rand_diag(::Type{T}; emin=-200, emax=200) where {T}
    x = rand(T)
    rand(Bool) || (x = -x)
    return T(ldexp(big(x), rand(emin:emax)))
end

function relerr(z, x, y)
    ref = big(x) * big(y)
    iszero(ref) && return zero(BigFloat)
    return abs(big(z) - ref) / abs(ref)
end

function diagnose33(label, generator; n=20_000)
    u = BigFloat(2)^(-53)
    bound = 44u^3 + 176u^4
    bound_fail = 0
    p_fail = 0
    norm_fail = 0
    finite_fail = 0
    max_ratio = zero(BigFloat)
    first_bound = nothing
    first_p = nothing

    for _ in 1:n
        x = generator(Float64x3)
        y = generator(Float64x3)
        (iszero(x) || iszero(y)) && continue
        z = tw_prod33_fast(x, y)
        if !isfinite(z)
            finite_fail += 1
            continue
        end
        e = relerr(z, x, y)
        ratio = e / bound
        max_ratio = max(max_ratio, ratio)
        if e > bound
            bound_fail += 1
            first_bound === nothing && (first_bound = (x._limbs, y._limbs, z._limbs, e, ratio))
        end
        if !paper_p_nonoverlap(z._limbs)
            p_fail += 1
            first_p === nothing && (first_p = (x._limbs, y._limbs, z._limbs))
        end
        MultiFloats.isnormalized(z) || (norm_fail += 1)
    end

    println("33/$label: bound_fail=$bound_fail p_fail=$p_fail norm_fail=$norm_fail finite_fail=$finite_fail max_err_over_bound=$(Float64(max_ratio))")
    first_bound === nothing || println("33/$label first_bound=", first_bound)
    first_p === nothing || println("33/$label first_p=", first_p)
end

function diagnose23(label, generator; n=20_000)
    u = BigFloat(2)^(-53)
    bound = 18u^3 + 75u^4
    bound_fail = 0
    p_fail = 0
    norm_fail = 0
    finite_fail = 0
    max_ratio = zero(BigFloat)
    first_bound = nothing
    first_p = nothing

    for _ in 1:n
        x = generator(Float64x2)
        y = generator(Float64x3)
        (iszero(x) || iszero(y)) && continue
        z = tw_prod23_fast(x, y)
        if !isfinite(z)
            finite_fail += 1
            continue
        end
        e = relerr(z, x, y)
        ratio = e / bound
        max_ratio = max(max_ratio, ratio)
        if e > bound
            bound_fail += 1
            first_bound === nothing && (first_bound = (x._limbs, y._limbs, z._limbs, e, ratio))
        end
        if !paper_p_nonoverlap(z._limbs)
            p_fail += 1
            first_p === nothing && (first_p = (x._limbs, y._limbs, z._limbs))
        end
        MultiFloats.isnormalized(z) || (norm_fail += 1)
    end

    println("23/$label: bound_fail=$bound_fail p_fail=$p_fail norm_fail=$norm_fail finite_fail=$finite_fail max_err_over_bound=$(Float64(max_ratio))")
    first_bound === nothing || println("23/$label first_bound=", first_bound)
    first_p === nothing || println("23/$label first_p=", first_p)
end

Random.seed!(0x3a19_d1a6)
setprecision(BigFloat, 1024) do
    ordinary(T) = rand(T) + T(0.5)
    wide(T) = wide_rand_diag(T)
    println("Triple-word 2019 non-failing diagnostics")
    println("CPU: ", Sys.CPU_NAME)
    diagnose33("ordinary", ordinary)
    diagnose33("wide", wide)
    diagnose23("ordinary", ordinary)
    diagnose23("wide", wide)
end
