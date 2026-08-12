# M6 correctness-first Float64x5 division candidate.
#
# Use the accepted x5 reciprocal only to form an initial quotient. Then perform
# one Karp-Markstein-style quotient correction using DIRECT x5 FMA for both the
# quotient residual and the correction update:
#
#     q0 = round_x5(x * inv(y))
#     r  = x - q0*y
#     q1 = q0 + r*inv(y)
#
# The direct residual/update avoid the rounded mul-then-add behavior rejected by
# M5. `reference_div` remains the independent acceptance authority.

export div5_safe

const _Float64x5_div5 = MultiFloat{Float64,5}

function _check_div5_input(x::_Float64x5_div5, y::_Float64x5_div5)
    isfinite(x) || throw(DomainError(x, "div5_safe requires a finite numerator"))
    isfinite(y) || throw(DomainError(y, "div5_safe requires a finite denominator"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "div5_safe requires a normalized numerator",
    ))
    MultiFloats.isnormalized(y) || throw(ArgumentError(
        "div5_safe requires a normalized denominator",
    ))
    iszero(y) && throw(DomainError(y, "div5_safe requires a nonzero denominator"))
    return nothing
end

function _div5_seed(x::_Float64x5_div5, y::_Float64x5_div5)
    _check_div5_input(x, y)
    iszero(x) && return zero(_Float64x5_div5), recip5_safe(y)
    invy = recip5_safe(y)
    q0 = mul5_safe(x, invy)
    (isfinite(q0) && MultiFloats.isnormalized(q0)) || throw(DomainError(
        (x, y),
        "div5_safe initial quotient is not finite and normalized",
    ))
    return q0, invy
end

function _div5_correct_once(
    x::_Float64x5_div5,
    y::_Float64x5_div5,
    q0::_Float64x5_div5,
    invy::_Float64x5_div5,
)
    # r = x - q0*y as one direct FMA.
    nq0 = _neg_recip5_safe(q0)
    residual = fma5_safe(nq0, y, x)

    # q1 = q0 + residual*invy as one direct FMA.
    q1 = fma5_safe(residual, invy, q0)
    (isfinite(q1) && MultiFloats.isnormalized(q1)) || throw(DomainError(
        (x, y),
        "div5_safe quotient correction produced an invalid expansion",
    ))
    return q1, residual
end

function _div5_uncorrected(x::_Float64x5_div5, y::_Float64x5_div5)
    q0, _ = _div5_seed(x, y)
    return q0
end

function _div5_one_correction(x::_Float64x5_div5, y::_Float64x5_div5)
    _check_div5_input(x, y)
    iszero(x) && return zero(_Float64x5_div5)
    q0, invy = _div5_seed(x, y)
    q1, _ = _div5_correct_once(x, y, q0, invy)
    return q1
end

"""
    div5_safe(x, y)

Correctness-first Float64x5 division candidate. It forms an initial quotient from
`recip5_safe(y)`, then performs one direct-FMA quotient-residual correction.

The current domain is the intersection of normalized finite x5 inputs,
`recip5_safe`'s denominator domain, and the conservative safe mul/FMA
product-underflow domains.
"""
function div5_safe(x::_Float64x5_div5, y::_Float64x5_div5)
    return _div5_one_correction(x, y)
end
