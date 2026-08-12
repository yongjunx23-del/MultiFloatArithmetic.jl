# Correctness-first x5-x8 arithmetic oracle.
#
# This file lives inside `MultiFloatArithmetic.Experimental`. It deliberately
# uses exact Rational{BigInt} arithmetic and then packs the exact dyadic result
# back into a normalized MultiFloat expansion. It is not a performance path.
# Optimized higher-limb networks must compare against this implementation before
# they can be considered for promotion.

@inline function _check_reference_format(::Type{T}, ::Val{N}) where {T,N}
    (T === Float32 || T === Float64) || throw(ArgumentError(
        "reference arithmetic currently supports Float32 and Float64 limbs",
    ))
    N in 5:8 || throw(ArgumentError(
        "reference arithmetic currently supports N = 5, 6, 7, or 8",
    ))
    return nothing
end

@inline function _reference_precision(::Type{T}, ::Val{N}) where {T,N}
    # `span` is the number of binary places from the largest finite exponent to
    # the smallest subnormal bit (inclusive). A product can approximately
    # double that span. The N-dependent and fixed margins cover carry growth and
    # the final rational-to-BigFloat conversion without relying on ambient MPFR
    # precision.
    span = exponent(floatmax(T)) - exponent(floatmin(T)) + precision(T)
    return 2 * span + 2 * ndigits(N; base=2) + 32
end

function _check_reference_input(x::MultiFloat{T,N}) where {T,N}
    _check_reference_format(T, Val{N}())
    isfinite(x) || throw(DomainError(x, "reference arithmetic requires finite inputs"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "reference arithmetic requires normalized MultiFloat inputs",
    ))
    return nothing
end

@inline _exact_value(x::MultiFloat) = Rational{BigInt}(x)

function _pack_exact(
    ::Type{MultiFloat{T,N}},
    value::Rational{BigInt},
) where {T,N}
    _check_reference_format(T, Val{N}())
    p = _reference_precision(T, Val{N}())
    return setprecision(BigFloat, p) do
        # All values entering this package are dyadic rationals. At the chosen
        # precision this conversion is exact for the x5-x8 Float32/Float64
        # domain; the test suite cross-checks it against a much larger MPFR
        # precision before this function is accepted as an oracle.
        MultiFloat{T,N}(BigFloat(value))
    end
end

"""
    reference_add(x, y)

Exact-rational reference addition for 5- through 8-limb Float32/Float64
`MultiFloat` values. This operation is intentionally allocation-heavy and must
not be used as a performance baseline.
"""
function reference_add(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    return _pack_exact(MultiFloat{T,N}, _exact_value(x) + _exact_value(y))
end

"""
    reference_sub(x, y)

Exact-rational reference subtraction for the x5-x8 research domain.
"""
function reference_sub(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    return _pack_exact(MultiFloat{T,N}, _exact_value(x) - _exact_value(y))
end

"""
    reference_mul(x, y)

Exact-rational reference multiplication for the x5-x8 research domain.
"""
function reference_mul(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    return _pack_exact(MultiFloat{T,N}, _exact_value(x) * _exact_value(y))
end

"""
    reference_fma(x, y, c)

Compute the exact dyadic value `x*y + c`, then pack it once into an N-limb
MultiFloat expansion. This is the primary correctness oracle for future x5-x8
direct FMA research.
"""
function reference_fma(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
    c::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    _check_reference_input(c)
    exact = _exact_value(x) * _exact_value(y) + _exact_value(c)
    return _pack_exact(MultiFloat{T,N}, exact)
end

# Lane-wise vector wrappers. These are intentionally defined through scalar
# extraction so the SIMD research code has an independent oracle rather than a
# second vectorized arithmetic network with shared failure modes.
function reference_add(
    x::MultiFloatVec{W,T,N},
    y::MultiFloatVec{W,T,N},
) where {W,T,N}
    scalars = ntuple(i -> reference_add(x[i], y[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end

function reference_sub(
    x::MultiFloatVec{W,T,N},
    y::MultiFloatVec{W,T,N},
) where {W,T,N}
    scalars = ntuple(i -> reference_sub(x[i], y[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end

function reference_mul(
    x::MultiFloatVec{W,T,N},
    y::MultiFloatVec{W,T,N},
) where {W,T,N}
    scalars = ntuple(i -> reference_mul(x[i], y[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end

function reference_fma(
    x::MultiFloatVec{W,T,N},
    y::MultiFloatVec{W,T,N},
    c::MultiFloatVec{W,T,N},
) where {W,T,N}
    scalars = ntuple(i -> reference_fma(x[i], y[i], c[i]), Val{W}())
    return MultiFloatVec{W,T,N}(scalars)
end
