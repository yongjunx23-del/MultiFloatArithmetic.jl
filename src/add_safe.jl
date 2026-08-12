# Correctness-first Float64x5-x8 addition baseline.
#
# Each same-order limb pair is summed with general TwoSum, producing an exact
# 2N-term expansion. The full expansion is renormalized before retaining and
# renormalizing its N-limb head. No FastTwoSum magnitude assumptions are used.

@inline function _add_safe_terms(
    x::NTuple{N,T},
    y::NTuple{N,T},
) where {N,T}
    pairs = ntuple(i -> two_sum(x[i], y[i]), Val(N))
    return ntuple(Val(2 * N)) do j
        i = (j + 1) >>> 1
        isodd(j) ? pairs[i][1] : pairs[i][2]
    end
end

@inline function _add_safe_full_limbs(
    x::NTuple{N,T},
    y::NTuple{N,T},
) where {N,T}
    return renormalize(_add_safe_terms(x, y))
end

@inline function _add_safe_tail_limbs(
    x::NTuple{N,T},
    y::NTuple{N,T},
) where {N,T}
    full = _add_safe_full_limbs(x, y)
    return ntuple(i -> full[i + N], Val(N))
end

"""
    add_safe_limbs(x, y)

Correctness-first Float64 x5-x8 addition on normalized limb tuples. The exact
2N-term TwoSum expansion is fully renormalized before its N-limb head is kept.

This is an iterative research baseline, not a fixed-cost performance kernel.
"""
@inline function add_safe_limbs(
    x::NTuple{N,Float64},
    y::NTuple{N,Float64},
) where {N}
    N in 5:8 || throw(ArgumentError(
        "add_safe_limbs currently supports N = 5, 6, 7, or 8",
    ))
    full = _add_safe_full_limbs(x, y)
    head = ntuple(i -> full[i], Val(N))
    return renormalize(head)
end

function _check_add_safe_input(x::MultiFloat{Float64,N}) where {N}
    N in 5:8 || throw(ArgumentError(
        "add_safe currently supports N = 5, 6, 7, or 8",
    ))
    isfinite(x) || throw(DomainError(x, "add_safe requires finite inputs"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "add_safe requires normalized MultiFloat inputs",
    ))
    return nothing
end

"""
    add_safe(x, y)

Safe scalar Float64 x5-x8 addition baseline. Future optimized addition networks
must match or improve this candidate's correctness/error contract before
performance comparison.
"""
function add_safe(
    x::MultiFloat{Float64,N},
    y::MultiFloat{Float64,N},
) where {N}
    _check_add_safe_input(x)
    _check_add_safe_input(y)
    return MultiFloat{Float64,N}(add_safe_limbs(x._limbs, y._limbs))
end

# Stable width-specific research names. These wrappers make width-scoped A/B
# studies explicit while sharing one correctness implementation.
@inline add5_safe(x::MultiFloat{Float64,5}, y::MultiFloat{Float64,5}) = add_safe(x, y)
@inline add6_safe(x::MultiFloat{Float64,6}, y::MultiFloat{Float64,6}) = add_safe(x, y)
@inline add7_safe(x::MultiFloat{Float64,7}, y::MultiFloat{Float64,7}) = add_safe(x, y)
@inline add8_safe(x::MultiFloat{Float64,8}, y::MultiFloat{Float64,8}) = add_safe(x, y)

@inline add5_safe_limbs(x::NTuple{5,Float64}, y::NTuple{5,Float64}) = add_safe_limbs(x, y)
@inline add6_safe_limbs(x::NTuple{6,Float64}, y::NTuple{6,Float64}) = add_safe_limbs(x, y)
@inline add7_safe_limbs(x::NTuple{7,Float64}, y::NTuple{7,Float64}) = add_safe_limbs(x, y)
@inline add8_safe_limbs(x::NTuple{8,Float64}, y::NTuple{8,Float64}) = add_safe_limbs(x, y)
