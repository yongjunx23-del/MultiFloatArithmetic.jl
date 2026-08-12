# M6 non-dyadic correctness oracles.
#
# Reciprocal and division differ from the M2 add/mul/FMA oracle: exact dyadic
# inputs generally produce a non-dyadic rational result. We therefore keep the
# source expression exact as Rational{BigInt}, round it at an explicit large MPFR
# precision, and only then pack once into MultiFloat. The permanent tests compare
# this guard precision against a much larger independent MPFR precision.
#
# This is deliberately allocation-heavy and is not a candidate arithmetic path.

@inline function _reference_nondyadic_precision(::Type{T}, ::Val{N}) where {T,N}
    _check_reference_format(T, Val{N}())
    # Much larger than the N-limb target precision (N*precision(T)). The exact
    # value can be arbitrarily close to a target rounding boundary, so this is a
    # high-guard empirical oracle rather than a proof of correct rounding. CI
    # cross-checks every accepted corpus against 16384-bit MPFR packing.
    return max(4096, 8 * N * precision(T) + 256)
end

function _pack_rounded_rational(
    ::Type{MultiFloat{T,N}},
    value::Rational{BigInt},
) where {T,N}
    _check_reference_format(T, Val{N}())
    p = _reference_nondyadic_precision(T, Val{N}())
    result = setprecision(BigFloat, p) do
        MultiFloat{T,N}(BigFloat(value))
    end
    isfinite(result) || throw(DomainError(
        value,
        "reference non-dyadic result overflowed the target MultiFloat format",
    ))
    if !iszero(value) && iszero(result)
        throw(DomainError(
            value,
            "reference non-dyadic result underflowed to zero in the target format",
        ))
    end
    MultiFloats.isnormalized(result) || error(
        "internal reference non-dyadic packing produced an unnormalized result",
    )
    return result
end

"""
    reference_recip(x)

High-guard MPFR reference reciprocal for normalized finite x5-x8 Float32/Float64
`MultiFloat` values. The exact input is first converted to `Rational{BigInt}`;
`1/x` remains exact as a rational until one explicit high-precision MPFR rounding
and final MultiFloat pack.

This is an Experimental correctness oracle, not a performance implementation.
"""
function reference_recip(x::MultiFloat{T,N}) where {T,N}
    _check_reference_input(x)
    qx = _exact_value(x)
    iszero(qx) && throw(DomainError(x, "reference_recip requires a nonzero input"))
    return _pack_rounded_rational(MultiFloat{T,N}, inv(qx))
end

"""
    reference_div(x, y)

High-guard MPFR reference division for normalized finite x5-x8 Float32/Float64
`MultiFloat` values. `x/y` is formed exactly as `Rational{BigInt}` before the
single high-precision MPFR rounding/pack.
"""
function reference_div(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    qy = _exact_value(y)
    iszero(qy) && throw(DomainError(y, "reference_div requires a nonzero divisor"))
    value = _exact_value(x) / qy
    # Exact zero numerator must remain exactly zero.
    iszero(value) && return zero(MultiFloat{T,N})
    return _pack_rounded_rational(MultiFloat{T,N}, value)
end

# Deliberately lane-wise vector reference wrappers, preserving implementation
# independence from future SIMD reciprocal/division candidates.
function reference_recip(x::MultiFloatVec{W,T,N}) where {W,T,N}
    scalars = ntuple(i -> reference_recip(x[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end

function reference_div(
    x::MultiFloatVec{W,T,N},
    y::MultiFloatVec{W,T,N},
) where {W,T,N}
    scalars = ntuple(i -> reference_div(x[i], y[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end
