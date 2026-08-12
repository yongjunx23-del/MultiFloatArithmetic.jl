"""
    MultiFloatArithmetic

Verification-oriented arithmetic kernels for fixed-length `MultiFloats.jl`
expansions. The top-level API contains only the empirically validated fused
multiply-add research kernel. Rejected, reference-only, or not-yet-accepted
candidates live in [`MultiFloatArithmetic.Experimental`](@ref).
"""
module MultiFloatArithmetic

using MultiFloats
import MultiFloats: MultiFloat, MultiFloatVec, fast_two_sum, two_prod, two_sum

export Experimental, MFLinearAlgebra, fma_fast, fma_fast_limbs

"""
    fma_fast_limbs(x, y, c)

Compute the fixed-length expansion candidate for `x*y + c` when `x`, `y`, and
`c` are normalized 2-, 3-, or 4-limb tuples.

This is an operand-relative research kernel. It is not correctly rounded and it
does not provide a strong result-relative guarantee under destructive
cancellation. Inputs are expected to be finite normalized expansions evaluated
under IEEE round-to-nearest semantics without algebraic reassociation.
"""

# Research kernels for fused multiply-add on fixed-length MultiFloat expansions.
# T may be Float32/Float64 or the SIMD lane type used by MultiFloatVec. Do not
# algebraically reorder these networks without re-verifying the resulting
# floating-point accumulation network.

@inline function fma_fast_limbs(
    x::NTuple{2,T},
    y::NTuple{2,T},
    c::NTuple{2,T},
) where {T}
    x0, x1 = x
    y0, y1 = y
    c0, c1 = c

    p00, e00 = two_prod(x0, y0)
    p01 = x0 * y1
    p10 = x1 * y0

    cross = p01 + p10
    low = e00 + c1
    low = low + cross

    high, carry = two_sum(p00, c0)
    carry = carry + low
    z0, z1 = fast_two_sum(high, carry)
    return (z0, z1)
end

@inline function fma_fast_limbs(
    x::NTuple{3,T},
    y::NTuple{3,T},
    c::NTuple{3,T},
) where {T}
    x0, x1, x2 = x
    y0, y1, y2 = y
    c0, c1, c2 = c

    p00, e00 = two_prod(x0, y0)
    p01, e01 = two_prod(x0, y1)
    p10, e10 = two_prod(x1, y0)

    p02 = x0 * y2
    p11 = x1 * y1
    p20 = x2 * y0

    sigma = (p02 + p20) + p11
    tail = ((e01 + e10) + sigma) + c2

    a, q1 = two_sum(p01, p10)
    a, q2 = two_sum(a, e00)
    a, q3 = two_sum(a, c1)
    tail = tail + ((q1 + q2) + q3)

    b, r = two_sum(p00, c0)
    m1, m2 = two_sum(r, a)
    m2 = m2 + tail

    w0, w1 = fast_two_sum(b, m1)
    w1, w2 = two_sum(w1, m2)
    z0, rho = two_sum(w0, w1)
    z1, z2 = fast_two_sum(rho, w2)

    # FPANVerifier refuted the leading non-overlap relation after the first
    # compression. A/B verification showed that one additional TwoSum followed
    # by a tail FastTwoSum is the smallest fixed-cost repair that proves both
    # output non-overlap relations. Do not remove or reorder these operations.
    z0, z1 = two_sum(z0, z1)
    z1, z2 = fast_two_sum(z1, z2)
    return (z0, z1, z2)
end

@inline function fma_fast_limbs(
    x::NTuple{4,T},
    y::NTuple{4,T},
    c::NTuple{4,T},
) where {T}
    x0, x1, x2, x3 = x
    y0, y1, y2, y3 = y
    c0, c1, c2, c3 = c

    p00, e00 = two_prod(x0, y0)
    p01, e01 = two_prod(x0, y1)
    p10, e10 = two_prod(x1, y0)
    p02, e02 = two_prod(x0, y2)
    p11, e11 = two_prod(x1, y1)
    p20, e20 = two_prod(x2, y0)

    p03 = x0 * y3
    p12 = x1 * y2
    p21 = x2 * y1
    p30 = x3 * y0
    diagonal3 = (p03 + p30) + (p12 + p21)

    b, r = two_sum(p00, c0)

    a1, f1 = two_sum(p01, p10)
    a1, f2 = two_sum(a1, e00)
    a1, f3 = two_sum(a1, c1)
    a1, f4 = two_sum(a1, r)

    a2, g1 = two_sum(p02, p20)
    a2, g2 = two_sum(a2, p11)
    e01e10, g4 = two_sum(e01, e10)
    a2, g3 = two_sum(a2, e01e10)
    a2, g5 = two_sum(a2, c2)
    a2, g6 = two_sum(a2, f1)
    a2, g7 = two_sum(a2, f2)
    a2, g8 = two_sum(a2, f3)
    a2, g9 = two_sum(a2, f4)

    t1 = e02 + e20
    t2 = e11 + diagonal3
    t3 = (t1 + t2) + c3

    t1 = g1 + g2
    t2 = g3 + g4
    t1 = t1 + t2
    t2 = g6 + g7
    t4 = g8 + g9
    t2 = t2 + t4
    t1 = t1 + t2
    t1 = t1 + g5
    a3 = t3 + t1

    # Current arXiv:2607.11391 QW normalization: five fixed cascades over the
    # four live words. FPANVerifier's FastTwoSum contract is exponent-order, not
    # the stronger textbook |a| >= |b| condition. The gate placement below is
    # the paper's proved fixed point: FTT / TTF / TFF / FFF / FFF. The former
    # two-pass 146-flop variant preserves the component sum but does not carry a
    # machine-proved full non-overlap guarantee under deep cancellation.

    # Pass 1: F T T
    w0, w1 = fast_two_sum(b, a1)
    w1, w2 = two_sum(w1, a2)
    w2, w3 = two_sum(w2, a3)

    # Pass 2: T T F
    w0, w1 = two_sum(w0, w1)
    w1, w2 = two_sum(w1, w2)
    w2, w3 = fast_two_sum(w2, w3)

    # Pass 3: T F F
    w0, w1 = two_sum(w0, w1)
    w1, w2 = fast_two_sum(w1, w2)
    w2, w3 = fast_two_sum(w2, w3)

    # Pass 4: F F F
    w0, w1 = fast_two_sum(w0, w1)
    w1, w2 = fast_two_sum(w1, w2)
    w2, w3 = fast_two_sum(w2, w3)

    # Pass 5: F F F
    z0, w1 = fast_two_sum(w0, w1)
    z1, w2 = fast_two_sum(w1, w2)
    z2, z3 = fast_two_sum(w2, w3)
    return (z0, z1, z2, z3)
end

"""
    fma_fast(x, y, c)

Evaluate the 2-, 3-, or 4-limb fused multiply-add research kernel for scalar
`MultiFloat` values or lane-wise for `MultiFloatVec` values. All three widths use
fixed-cost arithmetic networks; the x4 path uses the five-pass QW normalization
reproduced from the current arXiv:2607.11391 reference implementation.

See `docs/NUMERICAL_CONTRACT.md` before using this operation in residual,
refinement, stopping-criterion, or certificate code.
"""

# Keep explicit x4 entry points in addition to the generic width dispatch. The
# QW network is large enough that Julia 1.10 may otherwise fail to inline the
# generic `where N` wrapper into scalar array loops, inhibiting the same loop
# optimization obtained by the source-identical width-specialized audit kernel.
# This does not change arithmetic or the limb-level network.
@inline function fma_fast(
    x::MultiFloat{T,4},
    y::MultiFloat{T,4},
    c::MultiFloat{T,4},
) where {T}
    return MultiFloat{T,4}(fma_fast_limbs(x._limbs, y._limbs, c._limbs))
end

@inline function fma_fast(
    x::MultiFloatVec{W,T,4},
    y::MultiFloatVec{W,T,4},
    c::MultiFloatVec{W,T,4},
) where {W,T}
    return MultiFloatVec{W,T,4}(fma_fast_limbs(x._limbs, y._limbs, c._limbs))
end

@inline function fma_fast(
    x::MultiFloat{T,N},
    y::MultiFloat{T,N},
    c::MultiFloat{T,N},
) where {T,N}
    N in (2, 3, 4) || throw(ArgumentError("fma_fast currently supports N = 2, 3, or 4"))
    return MultiFloat{T,N}(fma_fast_limbs(x._limbs, y._limbs, c._limbs))
end

@inline function fma_fast(
    x::MultiFloatVec{W,T,N},
    y::MultiFloatVec{W,T,N},
    c::MultiFloatVec{W,T,N},
) where {W,T,N}
    N in (2, 3, 4) || throw(ArgumentError("fma_fast currently supports N = 2, 3, or 4"))
    return MultiFloatVec{W,T,N}(
        fma_fast_limbs(x._limbs, y._limbs, c._limbs),
    )
end

"""
    MultiFloatArithmetic.Experimental

Research candidates and correctness oracles that have not passed the acceptance
gates for the top-level API. Their names, behavior, and presence may change
without deprecation during the 0.x series.
"""
module Experimental

using MultiFloats
import MultiFloats: MultiFloat, MultiFloatVec, div_r, fast_two_sum, mfadd,
    renormalize, two_prod, two_sum

export add_safe, add_safe_limbs
export add5_safe, add5_safe_limbs, add6_safe, add6_safe_limbs
export add7_safe, add7_safe_limbs, add8_safe, add8_safe_limbs
export div_digits, div_digits_limbs
export fma5_safe, fma5_safe_limbs
export mul_safe, mul_safe_limbs
export mul5_safe, mul5_safe_limbs, mul6_safe, mul6_safe_limbs
export mul7_safe, mul7_safe_limbs, mul8_safe, mul8_safe_limbs
export mul_scalar, mul_scalar_limbs
export reference_add, reference_div, reference_fma, reference_mul, reference_recip,
    reference_sub

include("mul_scalar.jl")
include("division_digits.jl")
include("reference_arithmetic.jl")
include("reference_nondyadic.jl")
include("add_safe.jl")
include("mul_safe.jl")
include("fma5_safe.jl")

end # module Experimental

include("linear_algebra.jl")

end # module MultiFloatArithmetic
