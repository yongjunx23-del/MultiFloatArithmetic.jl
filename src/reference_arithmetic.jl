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
        MultiFloat{T,N}(BigFloat(value))
    end
end

function reference_add(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    return _pack_exact(MultiFloat{T,N}, _exact_value(x) + _exact_value(y))
end

function reference_sub(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    return _pack_exact(MultiFloat{T,N}, _exact_value(x) - _exact_value(y))
end

function reference_mul(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
) where {T,N}
    _check_reference_input(x)
    _check_reference_input(y)
    return _pack_exact(MultiFloat{T,N}, _exact_value(x) * _exact_value(y))
end

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

# M6 non-dyadic / irrational correctness oracles are loaded here so they share
# the exact input/domain helpers above but remain isolated from hot arithmetic.
include("reference_inverse_div_sqrt.jl")
