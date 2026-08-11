# Specialized expansion × one-limb multiplication.
#
# These routines are obtained by substituting y = (q, 0, ..., 0) into the
# corresponding MultiFloats v3 mfmul networks and deleting only product gates
# whose second operand is identically zero.  The accumulation/renormalization
# topology is intentionally retained in this first version.  This makes the
# candidate easy to validate directly against
#
#     mfmul(x, (q, 0, ..., 0), Val(N))
#
# before attempting more aggressive network simplification.

@inline function mul_scalar_limbs(
    x::NTuple{2,T},
    q::T,
    ::Val{2},
) where {T}
    p00, e00 = two_prod(x[1], q)
    p10 = x[2] * q
    e00 += p10
    p00, e00 = fast_two_sum(p00, e00)
    return (p00, e00)
end

@inline function mul_scalar_limbs(
    x::NTuple{3,T},
    q::T,
    ::Val{3},
) where {T}
    z = zero(T)

    p00, e00 = two_prod(x[1], q)
    p01 = z
    e01 = z
    p10, e10 = two_prod(x[2], q)
    p02 = z
    p11 = z
    p20 = x[3] * q

    # Preserve the mfmul(x, (q,0,0)) accumulation topology.
    p01, p10 = two_sum(p01, p10)
    e01 += e10
    p02 += p20
    e00, p01 = two_sum(e00, p01)
    p02 += p11
    p00, e00 = fast_two_sum(p00, e00)
    p01 += p10
    e01 += p02
    p01 += e01
    e00, p01 = two_sum(e00, p01)
    p00, e00 = fast_two_sum(p00, e00)
    e00, p01 = fast_two_sum(e00, p01)
    p00, e00 = fast_two_sum(p00, e00)
    return (p00, e00, p01)
end

@inline function mul_scalar_limbs(
    x::NTuple{4,T},
    q::T,
    ::Val{4},
) where {T}
    z = zero(T)

    p00, e00 = two_prod(x[1], q)
    p01 = z
    e01 = z
    p10, e10 = two_prod(x[2], q)
    p02 = z
    e02 = z
    p11 = z
    e11 = z
    p20, e20 = two_prod(x[3], q)
    p03 = z
    p12 = z
    p21 = z
    p30 = x[4] * q

    # Preserve the mfmul(x, (q,0,0,0)) accumulation topology.
    p01, p10 = two_sum(p01, p10)
    e01, e10 = two_sum(e01, e10)
    p02, p20 = two_sum(p02, p20)
    e02 += e20
    p03 += p30
    p12 += p21
    e00, p01 = two_sum(e00, p01)
    e01, p11 = two_sum(e01, p11)
    e10 += e02
    p20 += e11
    p03 += p12
    p00, e00 = fast_two_sum(p00, e00)
    p01, p10 = fast_two_sum(p01, p10)
    e01, p02 = two_sum(e01, p02)
    e10 += p03
    p11 += p20
    p01, e01 = two_sum(p01, e01)
    p10 += p11
    e10 += p02
    p10 += e01
    p01, p10 = two_sum(p01, p10)
    e00, p01 = two_sum(e00, p01)
    p10 += e10
    p00, e00 = fast_two_sum(p00, e00)
    p01, p10 = two_sum(p01, p10)
    e00, p01 = two_sum(e00, p01)
    p00, e00 = fast_two_sum(p00, e00)
    p01, p10 = fast_two_sum(p01, p10)
    e00, p01 = fast_two_sum(e00, p01)
    p00, e00 = fast_two_sum(p00, e00)
    p01, p10 = fast_two_sum(p01, p10)
    e00, p01 = fast_two_sum(e00, p01)
    p01, p10 = fast_two_sum(p01, p10)
    return (p00, e00, p01, p10)
end

@inline function mul_scalar(
    x::MultiFloat{T,N},
    q::T,
) where {T,N}
    N in (2, 3, 4) ||
        throw(ArgumentError("mul_scalar currently supports N = 2, 3, or 4"))
    return MultiFloat{T,N}(mul_scalar_limbs(x._limbs, q, Val{N}()))
end

@inline function mul_scalar(
    x::MultiFloatVec{W,T,N},
    q,
) where {W,T,N}
    N in (2, 3, 4) ||
        throw(ArgumentError("mul_scalar currently supports N = 2, 3, or 4"))
    return MultiFloatVec{W,T,N}(mul_scalar_limbs(x._limbs, q, Val{N}()))
end
