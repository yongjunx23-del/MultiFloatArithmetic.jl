# Correctness-first Float64x5-x8 multiplication family.
#
# For width N, every one of the N^2 limb products is evaluated with TwoProd so
# both rounded product and residual are retained. The resulting 2N^2 finite terms
# are canonicalized independently of operand order, fully renormalized, and only
# then truncated to N limbs.
#
# This is deliberately branchy/allocating and is not a performance path.

@inline _canonical_mul_term(x::Float64) = iszero(x) ? 0.0 : x

function _mul_safe_terms(
    x::NTuple{N,Float64},
    y::NTuple{N,Float64},
) where {N}
    N in 5:8 || throw(ArgumentError(
        "safe multiplication currently supports N = 5, 6, 7, or 8",
    ))
    terms = Vector{Float64}(undef, 2 * N^2)
    k = 1
    @inbounds for i in 1:N
        xi = x[i]
        for j in 1:N
            yj = y[j]
            p, e = two_prod(xi, yj)
            (isfinite(p) && isfinite(e)) || throw(DomainError(
                (xi, yj),
                "mul_safe requires finite TwoProd components",
            ))
            (iszero(p) || !issubnormal(p)) || throw(DomainError(
                (xi, yj),
                "mul_safe does not yet cover subnormal TwoProd products",
            ))
            (iszero(e) || !issubnormal(e)) || throw(DomainError(
                (xi, yj),
                "mul_safe does not yet cover subnormal TwoProd residuals",
            ))
            if iszero(p) && !iszero(xi) && !iszero(yj)
                throw(DomainError(
                    (xi, yj),
                    "mul_safe does not yet cover product underflow to zero",
                ))
            end
            terms[k] = _canonical_mul_term(p)
            terms[k + 1] = _canonical_mul_term(e)
            k += 2
        end
    end

    # The multiset of TwoProd components is invariant under x <-> y. Sorting by
    # a deterministic total key makes the subsequent iterative renormalization
    # bitwise commutative without relying on lexicographic pair order.
    sort!(terms; by=t -> (abs(t), t), rev=true)
    return Tuple(terms)
end

function _mul_safe_full_limbs(
    x::NTuple{N,Float64},
    y::NTuple{N,Float64},
) where {N}
    return renormalize(_mul_safe_terms(x, y))
end

# Backward-compatible internal names used by the already-frozen x5 diagnostic.
_mul5_terms(x::NTuple{5,Float64}, y::NTuple{5,Float64}) = _mul_safe_terms(x, y)
_mul5_full_limbs(x::NTuple{5,Float64}, y::NTuple{5,Float64}) =
    _mul_safe_full_limbs(x, y)

"""
    mul_safe_limbs(x, y)

Correctness-first Float64 multiplication baseline for 5- through 8-limb
expansions. It preserves all `2N^2` TwoProd product/residual components, fully
renormalizes the exact finite expansion, and returns a renormalized N-limb head.

The current contract excludes nonfinite and subnormal TwoProd components and
pair products that underflow to zero. This is not a fixed-cost performance
kernel.
"""
function mul_safe_limbs(
    x::NTuple{N,Float64},
    y::NTuple{N,Float64},
) where {N}
    N in 5:8 || throw(ArgumentError(
        "mul_safe_limbs currently supports N = 5, 6, 7, or 8",
    ))
    full = _mul_safe_full_limbs(x, y)
    head = ntuple(i -> full[i], Val{N}())
    return renormalize(head)
end

function _check_mul_safe_input(x::MultiFloat{Float64,N}) where {N}
    N in 5:8 || throw(ArgumentError(
        "mul_safe currently supports N = 5, 6, 7, or 8",
    ))
    isfinite(x) || throw(DomainError(x, "mul_safe requires finite inputs"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "mul_safe requires normalized MultiFloat inputs",
    ))
    return nothing
end

"""
    mul_safe(x, y)

Safe scalar Float64 multiplication research baseline for N=5:8. Future
higher-limb multiplication networks must match or improve its exact-oracle,
normalization, commutativity, and discarded-tail contract before performance
comparison.
"""
function mul_safe(
    x::MultiFloat{Float64,N},
    y::MultiFloat{Float64,N},
) where {N}
    _check_mul_safe_input(x)
    _check_mul_safe_input(y)
    return MultiFloat{Float64,N}(mul_safe_limbs(x._limbs, y._limbs))
end

mul5_safe_limbs(x::NTuple{5,Float64}, y::NTuple{5,Float64}) = mul_safe_limbs(x, y)
mul6_safe_limbs(x::NTuple{6,Float64}, y::NTuple{6,Float64}) = mul_safe_limbs(x, y)
mul7_safe_limbs(x::NTuple{7,Float64}, y::NTuple{7,Float64}) = mul_safe_limbs(x, y)
mul8_safe_limbs(x::NTuple{8,Float64}, y::NTuple{8,Float64}) = mul_safe_limbs(x, y)

mul5_safe(x::MultiFloat{Float64,5}, y::MultiFloat{Float64,5}) = mul_safe(x, y)
mul6_safe(x::MultiFloat{Float64,6}, y::MultiFloat{Float64,6}) = mul_safe(x, y)
mul7_safe(x::MultiFloat{Float64,7}, y::MultiFloat{Float64,7}) = mul_safe(x, y)
mul8_safe(x::MultiFloat{Float64,8}, y::MultiFloat{Float64,8}) = mul_safe(x, y)
