# Rejected default candidate: quotient-digit division.
#
# The outer algorithm extracts one base-limb quotient digit at a time. Residual
# products use the specialized expansion × one-limb network in mul_scalar.jl;
# residual subtraction still uses MultiFloats v3's mfadd network. The candidate
# is retained only to reproduce the correctness and performance study. Two
# hosted-runner A/Bs showed a consistently much slower scalar path and
# architecture-dependent Vec4 results. One Zen 3 x3 Vec4 case won, but the
# candidate did not establish a stable enough benefit to replace upstream `/`.

"""
    div_digits_limbs(x, y, Val(N))

Experimental quotient-digit division for normalized 2-, 3-, or 4-limb tuples.
The divisor must be finite and nonzero. This function is retained for research
reproducibility and is not an accepted replacement for upstream division.
"""
@inline function div_digits_limbs(
    x::NTuple{N,T},
    y::NTuple{N,T},
    ::Val{N},
) where {N,T}
    N in (2, 3, 4) ||
        throw(ArgumentError("div_digits currently supports N = 2, 3, or 4"))

    residual = x
    digits = ntuple(_ -> zero(T), Val{N}())

    @inbounds for k in 1:N
        qk = div_r(residual[1], y[1])
        digits = Base.setindex(digits, qk, k)

        product = mul_scalar_limbs(y, qk, Val{N}())
        residual = mfadd(residual, map(-, product), Val{N}())
    end

    return renormalize(digits)
end

"""
    div_digits(x, y)

Experimental scalar or lane-wise wrapper for [`div_digits_limbs`](@ref).
Use ordinary `x / y` for accepted package behavior.
"""
@inline function div_digits(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    return MultiFloat{T,N}(div_digits_limbs(x._limbs, y._limbs, Val{N}()))
end

@inline function div_digits(
    x::MultiFloatVec{W,T,N},
    y::MultiFloatVec{W,T,N},
) where {W,T,N}
    return MultiFloatVec{W,T,N}(
        div_digits_limbs(x._limbs, y._limbs, Val{N}()),
    )
end
