# Research candidate: quotient-digit division.
#
# This file deliberately changes only the OUTER division algorithm. Residual
# products/subtractions still use MultiFloats v3's existing mfmul/mfadd
# networks, so numerical behavior can be evaluated before introducing a new
# half-to-full multiplication primitive.

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

        # Reference-quality residual update. This zero-padded product currently
        # executes the full N×N mfmul network; a later research phase will
        # replace it with a verified N×1/half-to-full product and fused residual
        # subtraction if the quotient-digit outer algorithm proves worthwhile.
        qexp = _one_limb_expansion(qk, Val{N}())
        product = mfmul(y, qexp, Val{N}())
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
