# Correctness-first Float64x5 multiplication candidate.
#
# This is intentionally over-complete and branchy. Every one of the 25 limb
# products is evaluated with TwoProd so both the rounded product and residual are
# retained. The resulting 50 finite terms are canonicalized independently of
# operand order, fully renormalized, and only then truncated to five limbs.
#
# This is a rejection/reference baseline for future fixed-cost multiplication,
# not a hot-loop implementation.

@inline _canonical_mul_term(x::Float64) = iszero(x) ? 0.0 : x

function _mul5_terms(
    x::NTuple{5,Float64},
    y::NTuple{5,Float64},
)
    terms = Vector{Float64}(undef, 50)
    k = 1
    @inbounds for i in 1:5
        xi = x[i]
        for j in 1:5
            yj = y[j]
            p, e = two_prod(xi, yj)
            (isfinite(p) && isfinite(e)) || throw(DomainError(
                (xi, yj),
                "mul5_safe requires finite TwoProd components",
            ))
            # Keep the initial contract conservative around gradual underflow.
            # Zero residuals are ordinary; nonzero subnormal product/residual
            # components are rejected until a dedicated underflow proof exists.
            (iszero(p) || !issubnormal(p)) || throw(DomainError(
                (xi, yj),
                "mul5_safe does not yet cover subnormal TwoProd products",
            ))
            (iszero(e) || !issubnormal(e)) || throw(DomainError(
                (xi, yj),
                "mul5_safe does not yet cover subnormal TwoProd residuals",
            ))
            if iszero(p) && !iszero(xi) && !iszero(yj)
                throw(DomainError(
                    (xi, yj),
                    "mul5_safe does not yet cover product underflow to zero",
                ))
            end
            terms[k] = _canonical_mul_term(p)
            terms[k + 1] = _canonical_mul_term(e)
            k += 2
        end
    end

    # x*y and y*x produce the same multiset of TwoProd components. Canonical
    # sorting makes the subsequent iterative renormalization bitwise commutative
    # without assuming that lexicographic pair order is numerically harmless.
    sort!(terms; by=t -> (abs(t), t), rev=true)
    return Tuple(terms)
end

function _mul5_full_limbs(
    x::NTuple{5,Float64},
    y::NTuple{5,Float64},
)
    return renormalize(_mul5_terms(x, y))
end

"""
    mul5_safe_limbs(x, y)

Correctness-first Float64x5 multiplication baseline. It preserves every TwoProd
product/residual from the 5x5 limb-product grid, canonicalizes the 50-term exact
expansion, fully renormalizes it, then returns a renormalized five-limb head.

The current research contract excludes nonfinite TwoProd components and
subnormal/underflowing pair products. This is not a fixed-cost performance
kernel.
"""
function mul5_safe_limbs(
    x::NTuple{5,Float64},
    y::NTuple{5,Float64},
)
    full = _mul5_full_limbs(x, y)
    head = ntuple(i -> full[i], Val{5}())
    return renormalize(head)
end

function _check_mul5_input(x::MultiFloat{Float64,5})
    isfinite(x) || throw(DomainError(x, "mul5_safe requires finite inputs"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "mul5_safe requires normalized MultiFloat inputs",
    ))
    return nothing
end

"""
    mul5_safe(x, y)

Safe scalar Float64x5 multiplication research baseline. Future `mul5` networks
must match or improve its exact-oracle, normalization, commutativity, and
error-tail contract before performance comparison.
"""
function mul5_safe(
    x::MultiFloat{Float64,5},
    y::MultiFloat{Float64,5},
)
    _check_mul5_input(x)
    _check_mul5_input(y)
    return MultiFloat{Float64,5}(mul5_safe_limbs(x._limbs, y._limbs))
end
