# M6 correctness oracles for reciprocal, division, and square root.
#
# Unlike M2 add/mul/FMA, reciprocal/division results are generally non-dyadic
# and sqrt is generally irrational. Therefore a single finite BigFloat precision
# cannot be called exact. These references adaptively double MPFR precision until
# the final N-limb MultiFloat packing is stable at two consecutive precisions.
# CI independently repeats the pack at substantially larger precision.

export reference_inv, reference_div, reference_sqrt

@inline function _reference_start_precision(::Type{T}, ::Val{N}) where {T,N}
    _check_reference_format(T, Val{N}())
    return max(1024, 4 * N * precision(T) + 256)
end

function _stable_pack(
    ::Type{MultiFloat{T,N}},
    compute;
    max_precision::Int = 32768,
) where {T,N}
    p = _reference_start_precision(T, Val{N}())
    previous = nothing
    while p <= max_precision
        packed = setprecision(BigFloat, p) do
            value = compute()
            isfinite(value) || throw(DomainError(value, "reference result is not finite"))
            MultiFloat{T,N}(value)
        end
        MultiFloats.isnormalized(packed) || throw(ErrorException(
            "reference packing produced a non-normalized expansion at precision $p",
        ))
        if previous !== nothing && packed === previous
            return packed
        end
        previous = packed
        p *= 2
    end
    throw(ErrorException(
        "reference packing did not stabilize by $max_precision bits",
    ))
end

"""
    reference_inv(x)

Precision-stabilized reciprocal oracle for normalized finite Float32/Float64
5- through 8-limb values. Zero is rejected. The source value is the exact
`Rational{BigInt}` representation of `x`; only the non-dyadic reciprocal-to-
MultiFloat packing uses MPFR.
"""
function reference_inv(x::MultiFloat{T,N}) where {T,N}
    _check_reference_input(x)
    qx = _exact_value(x)
    iszero(qx) && throw(DomainError(x, "reference_inv requires a nonzero input"))
    return _stable_pack(MultiFloat{T,N}) do
        BigFloat(denominator(qx)) / BigFloat(numerator(qx))
    end
end

"""
    reference_div(x, y)

Precision-stabilized division oracle for normalized finite Float32/Float64
5- through 8-limb values. `y == 0` is rejected. The exact rational quotient is
formed conceptually before MPFR packing.
"""
function reference_div(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    qx = _exact_value(x)
    qy = _exact_value(y)
    iszero(qy) && throw(DomainError(y, "reference_div requires a nonzero denominator"))
    q = qx / qy
    return _stable_pack(MultiFloat{T,N}) do
        BigFloat(numerator(q)) / BigFloat(denominator(q))
    end
end

"""
    reference_sqrt(x)

Precision-stabilized square-root oracle for normalized finite Float32/Float64
5- through 8-limb values. Negative input is rejected; zero maps exactly to zero.
The generally irrational square root is evaluated independently at increasing
MPFR precision until its N-limb packing stabilizes.
"""
function reference_sqrt(x::MultiFloat{T,N}) where {T,N}
    _check_reference_input(x)
    qx = _exact_value(x)
    qx < 0 && throw(DomainError(x, "reference_sqrt requires a nonnegative input"))
    iszero(qx) && return zero(MultiFloat{T,N})
    return _stable_pack(MultiFloat{T,N}) do
        sqrt(BigFloat(numerator(qx)) / BigFloat(denominator(qx)))
    end
end

# Lane-wise wrappers remain intentionally scalar/oracle based so future vector
# Newton kernels cannot share their implementation failure modes.
function reference_inv(x::MultiFloatVec{W,T,N}) where {W,T,N}
    scalars = ntuple(i -> reference_inv(x[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end

function reference_div(
    x::MultiFloatVec{W,T,N},
    y::MultiFloatVec{W,T,N},
) where {W,T,N}
    scalars = ntuple(i -> reference_div(x[i], y[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end

function reference_sqrt(x::MultiFloatVec{W,T,N}) where {W,T,N}
    scalars = ntuple(i -> reference_sqrt(x[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end
