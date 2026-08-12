# Correctness-first Float64x5 addition candidate.
#
# This is intentionally over-complete and iterative. It avoids FastTwoSum
# magnitude assumptions entirely: pairwise TwoSum captures each same-order limb
# sum exactly, the resulting 10-term expansion is fully renormalized, and only
# then is the five-limb head retained. The discarded normalized tail is available
# to research diagnostics so an explicit error bound can be derived before any
# fixed-cost minimization.

@inline function _add5_terms(
    x::NTuple{5,T},
    y::NTuple{5,T},
) where {T}
    s1, e1 = two_sum(x[1], y[1])
    s2, e2 = two_sum(x[2], y[2])
    s3, e3 = two_sum(x[3], y[3])
    s4, e4 = two_sum(x[4], y[4])
    s5, e5 = two_sum(x[5], y[5])
    return (s1, e1, s2, e2, s3, e3, s4, e4, s5, e5)
end

@inline function _add5_full_limbs(
    x::NTuple{5,T},
    y::NTuple{5,T},
) where {T}
    return renormalize(_add5_terms(x, y))
end

"""
    add5_safe_limbs(x, y)

Correctness-first Float64 five-limb addition candidate. The inputs must be
normalized finite Float64 expansion limbs. The algorithm fully renormalizes an
exact 10-term TwoSum expansion before discarding the low five terms.

This is a research baseline, not a fixed-cost performance kernel.
"""
@inline function add5_safe_limbs(
    x::NTuple{5,Float64},
    y::NTuple{5,Float64},
)
    full = _add5_full_limbs(x, y)
    head = ntuple(i -> full[i], Val{5}())
    return renormalize(head)
end

@inline function _add5_tail_limbs(
    x::NTuple{5,T},
    y::NTuple{5,T},
) where {T}
    full = _add5_full_limbs(x, y)
    return ntuple(i -> full[i + 5], Val{5}())
end

function _check_add5_input(x::MultiFloat{Float64,5})
    isfinite(x) || throw(DomainError(x, "add5_safe requires finite inputs"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "add5_safe requires normalized MultiFloat inputs",
    ))
    return nothing
end

"""
    add5_safe(x, y)

Safe scalar Float64x5 addition research baseline. Future `add5` networks must
match or improve this candidate's correctness/error contract before performance
comparison.
"""
function add5_safe(
    x::MultiFloat{Float64,5},
    y::MultiFloat{Float64,5},
)
    _check_add5_input(x)
    _check_add5_input(y)
    return MultiFloat{Float64,5}(add5_safe_limbs(x._limbs, y._limbs))
end
