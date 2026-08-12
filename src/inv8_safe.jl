# M6 correctness-first Float64x8 reciprocal candidate.
#
# Use trusted upstream Float64x4 inverse as a ~212-bit seed, lift it to x8, then
# apply direct x8 Newton corrections. The update uses direct `fma_safe` rather
# than rounded `mul_safe` followed by `add_safe`; M5 demonstrated that such
# rounded composition can lose the direct result under cancellation.

export inv8_safe

const _Float64x4_inv8 = MultiFloat{Float64,4}
const _Float64x8_inv8 = MultiFloat{Float64,8}

@inline _canonical_neg_inv8_limb(x::Float64) = iszero(x) ? 0.0 : -x

function _neg_inv8_safe(x::_Float64x8_inv8)
    limbs = ntuple(i -> _canonical_neg_inv8_limb(x._limbs[i]), Val{8}())
    return _Float64x8_inv8(renormalize(limbs))
end

function _check_inv8_input(x::_Float64x8_inv8)
    isfinite(x) || throw(DomainError(x, "inv8_safe requires a finite input"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "inv8_safe requires a normalized MultiFloat input",
    ))
    iszero(x) && throw(DomainError(x, "inv8_safe requires a nonzero input"))

    lead = x._limbs[1]
    iszero(lead) && throw(DomainError(x, "inv8_safe requires a nonzero leading limb"))
    issubnormal(lead) && throw(DomainError(
        x,
        "inv8_safe does not yet cover subnormal leading limbs",
    ))
    lead_inv = inv(lead)
    (isfinite(lead_inv) && !iszero(lead_inv) && !issubnormal(lead_inv)) || throw(DomainError(
        x,
        "inv8_safe leading-limb reciprocal is outside the current normal seed domain",
    ))
    return nothing
end

function _inv8_seed_x4(x::_Float64x8_inv8)
    _check_inv8_input(x)
    x4 = _Float64x4_inv8(ntuple(i -> x._limbs[i], Val{4}()))
    MultiFloats.isnormalized(x4) || throw(ArgumentError(
        "inv8_safe requires the leading x4 truncation to remain normalized",
    ))
    y4 = inv(x4)
    (isfinite(y4) && MultiFloats.isnormalized(y4)) || throw(DomainError(
        x,
        "inv8_safe x4 reciprocal seed is not finite and normalized",
    ))
    limbs = ntuple(i -> i <= 4 ? y4._limbs[i] : 0.0, Val{8}())
    y8 = _Float64x8_inv8(renormalize(limbs))
    MultiFloats.isnormalized(y8) || error("internal x4-to-x8 reciprocal lift is not normalized")
    return y8
end

function _inv8_newton_once(x::_Float64x8_inv8, y::_Float64x8_inv8)
    minus_one = _Float64x8_inv8(-1.0)
    residual = fma_safe(x, y, minus_one)        # x*y - 1
    ny = _neg_inv8_safe(y)
    next_y = fma_safe(ny, residual, y)          # y - y*residual
    (isfinite(next_y) && MultiFloats.isnormalized(next_y)) || throw(DomainError(
        (x, y),
        "inv8_safe Newton correction produced an invalid expansion",
    ))
    return next_y, residual
end

function _inv8_direct_one_correction(x::_Float64x8_inv8)
    y0 = _inv8_seed_x4(x)
    y1, _ = _inv8_newton_once(x, y0)
    return y1
end

function _inv8_direct_two_corrections(x::_Float64x8_inv8)
    y0 = _inv8_seed_x4(x)
    y1, _ = _inv8_newton_once(x, y0)
    y2, _ = _inv8_newton_once(x, y1)
    return y2
end

# Diagnostic only until the A/B proves it is required. If y3 is bitwise equal to
# y2 in cases where y2 still differs from reference_inv, target-width Newton has
# reached a fixed point and more same-width corrections cannot select the correct
# neighboring x8 rounding.
function _inv8_direct_three_corrections(x::_Float64x8_inv8)
    y0 = _inv8_seed_x4(x)
    y1, _ = _inv8_newton_once(x, y0)
    y2, _ = _inv8_newton_once(x, y1)
    y3, _ = _inv8_newton_once(x, y2)
    return y3
end

"""
    inv8_safe(x)

Correctness-first Float64x8 reciprocal experiment. The public candidate remains
the two-correction path while CI compares one, two, and three direct-FMA Newton
corrections against the adaptive `reference_inv` oracle. It must not be promoted
until bitwise oracle agreement is established.
"""
function inv8_safe(x::_Float64x8_inv8)
    return _inv8_direct_two_corrections(x)
end
