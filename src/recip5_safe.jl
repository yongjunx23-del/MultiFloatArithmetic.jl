# M6 correctness-first Float64x5 reciprocal candidate.
#
# The accepted candidate uses a trusted upstream Float64x4 reciprocal as a
# ~212-bit seed, lifts that normalized value exactly to x5, then performs ONE x5
# Newton correction with direct safe FMA. The previous Float64-seed / three-x5-
# correction construction is retained only as an A/B diagnostic.
#
# If z approximates 1/a, define e = a*z - 1. Then
#
#     z_new = z - z*e = z*(2 - a*z).
#
# Both e and z_new are evaluated with direct fma5_safe, avoiding the rounded
# mul-then-add behavior M5 rejected under cancellation.

const _Float64x4_recip5 = MultiFloat{Float64,4}
const _Float64x5_recip5 = MultiFloat{Float64,5}

@inline function _canonical_neg_recip5_limb(x::Float64)
    return iszero(x) ? 0.0 : -x
end

function _neg_recip5_safe(x::_Float64x5_recip5)
    limbs = ntuple(i -> _canonical_neg_recip5_limb(x._limbs[i]), Val{5}())
    return _Float64x5_recip5(renormalize(limbs))
end

function _check_recip5_input(x::_Float64x5_recip5)
    isfinite(x) || throw(DomainError(x, "recip5_safe requires a finite input"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "recip5_safe requires a normalized MultiFloat input",
    ))
    iszero(x) && throw(DomainError(x, "recip5_safe requires a nonzero input"))

    lead = x._limbs[1]
    iszero(lead) && throw(DomainError(x, "recip5_safe requires a nonzero leading limb"))
    issubnormal(lead) && throw(DomainError(
        x,
        "recip5_safe does not yet cover subnormal leading limbs",
    ))
    lead_inv = inv(lead)
    (isfinite(lead_inv) && !iszero(lead_inv) && !issubnormal(lead_inv)) || throw(DomainError(
        x,
        "recip5_safe leading-limb reciprocal is outside the current normal seed domain",
    ))
    return nothing
end

function _recip5_seed_x4(x::_Float64x5_recip5)
    _check_recip5_input(x)
    x4 = _Float64x4_recip5(ntuple(i -> x._limbs[i], Val{4}()))
    MultiFloats.isnormalized(x4) || throw(ArgumentError(
        "recip5_safe requires the leading x4 truncation to remain normalized",
    ))

    # MultiFloats.jl's x4 inverse follows its Karp-Markstein precision-doubling
    # path. Treat it as the trusted lower-width seed; the x5 oracle remains the
    # independent acceptance authority.
    z4 = inv(x4)
    (isfinite(z4) && MultiFloats.isnormalized(z4)) || throw(DomainError(
        x,
        "recip5_safe x4 reciprocal seed is not finite and normalized",
    ))
    limbs = ntuple(i -> i <= 4 ? z4._limbs[i] : 0.0, Val{5}())
    z5 = _Float64x5_recip5(renormalize(limbs))
    MultiFloats.isnormalized(z5) || error("internal x4-to-x5 reciprocal lift is not normalized")
    return z5
end

function _recip5_seed_float64(x::_Float64x5_recip5)
    _check_recip5_input(x)
    return _Float64x5_recip5(inv(x._limbs[1]))
end

function _recip5_newton_once(x::_Float64x5_recip5, z::_Float64x5_recip5)
    minus_one = _Float64x5_recip5(-1.0)
    residual = fma5_safe(x, z, minus_one)        # x*z - 1
    nz = _neg_recip5_safe(z)
    next_z = fma5_safe(nz, residual, z)          # z - z*residual
    (isfinite(next_z) && MultiFloats.isnormalized(next_z)) || throw(DomainError(
        (x, z),
        "recip5_safe Newton correction produced an invalid expansion",
    ))
    return next_z, residual
end

# A/B diagnostic: Float64 seed followed by three x5 Newton corrections.
function _recip5_float64_three_corrections(x::_Float64x5_recip5)
    z = _recip5_seed_float64(x)
    for _ in 1:3
        z, _ = _recip5_newton_once(x, z)
    end
    return z
end

# A/B diagnostic and accepted candidate: x4 seed followed by one x5 correction.
function _recip5_x4_one_correction(x::_Float64x5_recip5)
    z0 = _recip5_seed_x4(x)
    z1, _ = _recip5_newton_once(x, z0)
    return z1
end

"""
    recip5_safe(x)

Correctness-first Float64x5 reciprocal candidate. It uses the trusted upstream
Float64x4 reciprocal as a precision-doubled seed and performs one direct-FMA x5
Newton correction.

The current domain requires a finite normalized nonzero x5 input whose leading
limb and leading-limb reciprocal are normal Float64 values. It also inherits the
safe direct-FMA product-side underflow exclusions.
"""
function recip5_safe(x::_Float64x5_recip5)
    return _recip5_x4_one_correction(x)
end
