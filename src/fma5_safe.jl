# Correctness-first Float64x5 direct FMA candidate.
#
# This combines the exact TwoProd product components with the five exact addend
# limbs before any N-limb truncation. It deliberately does not implement FMA as
# `mul_safe` followed by `add_safe`, so the direct-FMA rounding path has its own
# exact full-expansion accounting and can later be minimized independently.

@inline _canonical_fma_term(x::Float64) = iszero(x) ? 0.0 : x

function _fma5_terms(
    x::NTuple{5,Float64},
    y::NTuple{5,Float64},
    c::NTuple{5,Float64},
)
    # Reuse only the exact 25 TwoProd product/residual decomposition from M4,
    # not its five-limb multiplication result.
    product_terms = _mul_safe_terms(x, y)
    terms = Vector{Float64}(undef, 55)
    @inbounds for i in 1:50
        terms[i] = product_terms[i]
    end
    @inbounds for i in 1:5
        ci = c[i]
        isfinite(ci) || throw(DomainError(ci, "fma5_safe requires finite addend limbs"))
        terms[50 + i] = _canonical_fma_term(ci)
    end

    # Product-component multiset is invariant under x <-> y; c is unchanged.
    # Canonical sorting therefore gives bitwise x/y symmetry before iterative
    # renormalization.
    sort!(terms; by=t -> (abs(t), t), rev=true)
    return Tuple(terms)
end

function _fma5_full_limbs(
    x::NTuple{5,Float64},
    y::NTuple{5,Float64},
    c::NTuple{5,Float64},
)
    full = renormalize(_fma5_terms(x, y, c))
    all(isfinite, full) || throw(DomainError(
        (x, y, c),
        "fma5_safe full expansion overflowed during renormalization",
    ))
    return full
end

"""
    fma5_safe_limbs(x, y, c)

Correctness-first Float64x5 direct fused-multiply-add baseline. It forms a
55-component exact `x*y+c` expansion from all 25 TwoProd product/residual pairs
plus the five addend limbs, fully renormalizes it, then returns a renormalized
five-limb head.

This is not a fixed-cost performance kernel. Its product-side domain inherits the
current `mul_safe` exclusion of subnormal/underflowing TwoProd components.
"""
function fma5_safe_limbs(
    x::NTuple{5,Float64},
    y::NTuple{5,Float64},
    c::NTuple{5,Float64},
)
    full = _fma5_full_limbs(x, y, c)
    head = ntuple(i -> full[i], Val{5}())
    return renormalize(head)
end

function _check_fma5_input(x::MultiFloat{Float64,5})
    isfinite(x) || throw(DomainError(x, "fma5_safe requires finite inputs"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "fma5_safe requires normalized MultiFloat inputs",
    ))
    return nothing
end

"""
    fma5_safe(x, y, c)

Safe scalar Float64x5 direct-FMA research baseline. Future direct x5 fused
networks must match or improve its exact-oracle, normalization, x/y-symmetry,
and discarded-tail contract before performance comparison.
"""
function fma5_safe(
    x::MultiFloat{Float64,5},
    y::MultiFloat{Float64,5},
    c::MultiFloat{Float64,5},
)
    _check_fma5_input(x)
    _check_fma5_input(y)
    _check_fma5_input(c)
    return MultiFloat{Float64,5}(fma5_safe_limbs(x._limbs, y._limbs, c._limbs))
end

# Load the generic x5-x8 direct-FMA family after freezing the x5 wrapper above.
include("fma_safe.jl")

# M6 reciprocal and division use direct x5 FMA for correction steps.
include("recip5_safe.jl")
include("div5_safe.jl")
