# Correctness-first Float64x5-x8 direct FMA family.
#
# Width N combines every exact TwoProd product/residual component from the M4
# multiplication decomposition with all N addend limbs before any N-limb
# truncation. The full `2N^2 + N` component expansion is canonicalized, fully
# renormalized, and only then reduced to an N-limb head.
#
# This file is included after the already-frozen x5 implementation. It adds the
# generic family and x6-x8 wrappers without overriding the existing x5 wrapper;
# tests require generic x5 to be bitwise identical to `fma5_safe`.

export fma_safe, fma_safe_limbs
export fma6_safe, fma6_safe_limbs, fma7_safe, fma7_safe_limbs
export fma8_safe, fma8_safe_limbs

@inline _canonical_fma_safe_term(x::Float64) = iszero(x) ? 0.0 : x

function _fma_safe_terms(
    x::NTuple{N,Float64},
    y::NTuple{N,Float64},
    c::NTuple{N,Float64},
) where {N}
    N in 5:8 || throw(ArgumentError(
        "safe direct FMA currently supports N = 5, 6, 7, or 8",
    ))
    product_terms = _mul_safe_terms(x, y)
    nproduct = 2 * N^2
    terms = Vector{Float64}(undef, nproduct + N)
    @inbounds for i in 1:nproduct
        terms[i] = product_terms[i]
    end
    @inbounds for i in 1:N
        ci = c[i]
        isfinite(ci) || throw(DomainError(ci, "fma_safe requires finite addend limbs"))
        terms[nproduct + i] = _canonical_fma_safe_term(ci)
    end
    sort!(terms; by=t -> (abs(t), t), rev=true)
    return Tuple(terms)
end

function _fma_safe_full_limbs(
    x::NTuple{N,Float64},
    y::NTuple{N,Float64},
    c::NTuple{N,Float64},
) where {N}
    full = renormalize(_fma_safe_terms(x, y, c))
    all(isfinite, full) || throw(DomainError(
        (x, y, c),
        "fma_safe full expansion overflowed during renormalization",
    ))
    return full
end

"""
    fma_safe_limbs(x, y, c)

Correctness-first direct Float64 FMA baseline for 5- through 8-limb expansions.
It forms one exact `x*y+c` expansion from all `2N^2` TwoProd product/residual
components plus N addend limbs, fully renormalizes it, and returns a renormalized
N-limb head.

The product-side domain inherits `mul_safe`'s current exclusion of nonfinite,
subnormal, and underflow-to-zero TwoProd components. This is not a fixed-cost
performance kernel.
"""
function fma_safe_limbs(
    x::NTuple{N,Float64},
    y::NTuple{N,Float64},
    c::NTuple{N,Float64},
) where {N}
    N in 5:8 || throw(ArgumentError(
        "fma_safe_limbs currently supports N = 5, 6, 7, or 8",
    ))
    full = _fma_safe_full_limbs(x, y, c)
    head = ntuple(i -> full[i], Val{N}())
    return renormalize(head)
end

function _check_fma_safe_input(x::MultiFloat{Float64,N}) where {N}
    N in 5:8 || throw(ArgumentError(
        "fma_safe currently supports N = 5, 6, 7, or 8",
    ))
    isfinite(x) || throw(DomainError(x, "fma_safe requires finite inputs"))
    MultiFloats.isnormalized(x) || throw(ArgumentError(
        "fma_safe requires normalized MultiFloat inputs",
    ))
    return nothing
end

"""
    fma_safe(x, y, c)

Safe scalar direct Float64 FMA research baseline for N=5:8. Future fused
higher-limb networks must match or improve its exact-oracle, normalization,
x/y-symmetry, destructive-cancellation, and discarded-tail contract before
performance comparison.
"""
function fma_safe(
    x::MultiFloat{Float64,N},
    y::MultiFloat{Float64,N},
    c::MultiFloat{Float64,N},
) where {N}
    _check_fma_safe_input(x)
    _check_fma_safe_input(y)
    _check_fma_safe_input(c)
    return MultiFloat{Float64,N}(fma_safe_limbs(x._limbs, y._limbs, c._limbs))
end

fma6_safe_limbs(x::NTuple{6,Float64}, y::NTuple{6,Float64}, c::NTuple{6,Float64}) =
    fma_safe_limbs(x, y, c)
fma7_safe_limbs(x::NTuple{7,Float64}, y::NTuple{7,Float64}, c::NTuple{7,Float64}) =
    fma_safe_limbs(x, y, c)
fma8_safe_limbs(x::NTuple{8,Float64}, y::NTuple{8,Float64}, c::NTuple{8,Float64}) =
    fma_safe_limbs(x, y, c)

fma6_safe(x::MultiFloat{Float64,6}, y::MultiFloat{Float64,6}, c::MultiFloat{Float64,6}) =
    fma_safe(x, y, c)
fma7_safe(x::MultiFloat{Float64,7}, y::MultiFloat{Float64,7}, c::MultiFloat{Float64,7}) =
    fma_safe(x, y, c)
fma8_safe(x::MultiFloat{Float64,8}, y::MultiFloat{Float64,8}, c::MultiFloat{Float64,8}) =
    fma_safe(x, y, c)

# M6 reciprocal candidate depends on the accepted x8 direct-FMA family above.
include("inv8_safe.jl")
