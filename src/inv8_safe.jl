# Correctness-first Float64x8 reciprocal candidate.
#
# The first four limbs are inverted with trusted upstream Float64x4 arithmetic,
# whose Karp-Markstein implementation internally refines 1 -> 2 -> 4. That x4
# approximation is lifted exactly to x8 and refined with the accepted safe x8
# direct-FMA/multiplication/addition primitives. No BigFloat appears in the
# candidate iteration.

export inv8_safe

const _Float64x4 = MultiFloat{Float64,4}
const _Float64x8 = MultiFloat{Float64,8}

@inline function _canonical_neg_limb(x::Float64)
    return iszero(x) ? 0.0 : -x
end

function _neg_safe(x::MultiFloat{Float64,N}) where {N}
    limbs = ntuple(i -> _canonical_neg_limb(x._limbs[i]), Val{N}())
    return MultiFloat{Float64,N}(renormalize(limbs))
end

function _inv8_seed(x::_Float64x8)
    x4 = _Float64x4(ntuple(i -> x._limbs[i], Val{4}()))
    MultiFloats.isnormalized(x4) || throw(ArgumentError(
        "inv8_safe requires the leading x4 truncation to remain normalized",
    ))

    # MultiFloats.jl's x4 Base.inv path uses its Karp-Markstein mfinv recursion,
    # seeded from Float64 and refined through x2 to x4.
    y4 = inv(x4)
    (isfinite(y4) && MultiFloats.isnormalized(y4)) || throw(DomainError(
        x,
        "inv8_safe x4 reciprocal seed is not a finite normalized expansion",
    ))

    limbs = ntuple(i -> i <= 4 ? y4._limbs[i] : 0.0, Val{8}())
    return _Float64x8(renormalize(limbs))
end

function _inv8_newton_once(x::_Float64x8, y::_Float64x8)
    # Newton: y <- y - y*(x*y - 1).
    # fma_safe forms x*y-1 directly, avoiding an intermediate rounded product.
    residual = fma_safe(x, y, _Float64x8(-1.0))
    correction = mul_safe(y, residual)
    next_y = add_safe(y, _neg_safe(correction))
    (isfinite(next_y) && MultiFloats.isnormalized(next_y)) || throw(DomainError(
        (x, y),
        "inv8_safe Newton correction produced an invalid expansion",
    ))
    return next_y
end

function _check_inv8_input(x::_Float64x8)
    isfinite(x) || throw(DomainError(x, "inv8_safe requires a finite input"))
    iszero(x) && throw(DomainError(x, "inv8_safe requires a nonzero input"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "inv8_safe requires a normalized Float64x8 input",
    ))
    return nothing
end

"""
    inv8_safe(x)

Correctness-first Float64x8 reciprocal candidate. It uses the trusted upstream
x4 Karp-Markstein reciprocal as a 1->2->4 seed, lifts it to x8, then performs two
x8 Newton corrections using `fma_safe`, `mul_safe`, and `add_safe`.

The two-correction form is the initial acceptance baseline. Diagnostics compare
one versus two x8 corrections against `reference_inv` before any attempt to
remove the second correction.
"""
function inv8_safe(x::_Float64x8)
    _check_inv8_input(x)
    y0 = _inv8_seed(x)
    y1 = _inv8_newton_once(x, y0)
    return _inv8_newton_once(x, y1)
end

# Diagnostic entry point: one x8 correction after the same trusted x4 seed.
function _inv8_safe_one_correction(x::_Float64x8)
    _check_inv8_input(x)
    return _inv8_newton_once(x, _inv8_seed(x))
end
