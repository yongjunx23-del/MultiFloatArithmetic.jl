# Research candidate: quotient-digit division.
#
# The outer algorithm extracts one base-limb quotient digit at a time. Residual
# products use the specialized expansion × one-limb network in mul_scalar.jl;
# residual subtraction itself still uses MultiFloats v3's mfadd network. This
# keeps the strong-cancellation part on an existing upstream arithmetic path
# while removing the wasteful zero-padded N×N multiplication.

@inline _one_limb_expansion(q::T, ::Val{N}) where {T,N} =
    ntuple(i -> isone(i) ? q : zero(T), Val{N}())

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
