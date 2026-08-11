# Triple-word arithmetic research kernels from:
#
# N. Fabiano, J.-M. Muller, J. Picot,
# "Algorithms for Triple-Word Arithmetic",
# IEEE Transactions on Computers 68(11), 2019,
# DOI: 10.1109/TC.2019.2918451.
#
# This first port implements the fast multiplication primitives that are fully
# specified in the body of the paper:
#   * Algorithm 10: 3Prod_fast(3,3), relative error <= 44u^3 + 176u^4
#   * Algorithm 12: 3Prod_fast(2,3), relative error <= 18u^3 + 75u^4
#
# It is intentionally scalar-only. The paper's VSEB normalization contains a
# value-dependent test, so a lane-wise SIMD port needs a separate masked or
# branch-free design rather than blindly applying scalar control flow to
# MultiFloatVec.

# VecSum (Algorithm 4) specialized to the small arities used below.  The paper
# proves that Fast2Sum can replace 2Sum for these intermediate sequences.
@inline function _tw_vecsum3(a::Float64, b::Float64, c::Float64)
    e2 = c
    e1, e2 = fast_two_sum(b, e2)
    e0, e1 = fast_two_sum(a, e1)
    return e0, e1, e2
end

@inline function _tw_vecsum4(a::Float64, b::Float64, c::Float64, d::Float64)
    e3 = d
    e2, e3 = fast_two_sum(c, e3)
    e1, e2 = fast_two_sum(b, e2)
    e0, e1 = fast_two_sum(a, e1)
    return e0, e1, e2, e3
end

# VSEB(2) (Algorithm 5) for three F-nonoverlapping inputs.  Only the first two
# P-nonoverlapping outputs are retained, exactly as Algorithms 10 and 12 do.
@inline function _tw_vseb2_3(a::Float64, b::Float64, c::Float64)
    q, err = fast_two_sum(a, b)
    if !iszero(err)
        r0 = q
        q = err
        r1, _ = fast_two_sum(q, c)
        return r0, r1
    end
    r0, r1 = fast_two_sum(q, c)
    return r0, r1
end

"""
    tw_prod23_fast_limbs(x::NTuple{2,Float64}, y::NTuple{3,Float64})

Fabiano–Muller–Picot Algorithm 12, fast double-word × triple-word product.
The proven relative error bound (assuming the paper's TW/DW input contracts and
no under/overflow) is `18u^3 + 75u^4`, where `u = 2^-53` for binary64.
"""
@inline function tw_prod23_fast_limbs(
    x::NTuple{2,Float64},
    y::NTuple{3,Float64},
)
    x0, x1 = x
    y0, y1, y2 = y

    z00p, z00m = two_prod(x0, y0)
    z01p, z01m = two_prod(x0, y1)
    z10p, z10m = two_prod(x1, y0)

    b0, b1, b2 = _tw_vecsum3(z00m, z01p, z10p)
    c = fma(x1, y1, b2)

    z31 = fma(x0, y2, z10m)
    z3 = z31 + z01m
    s3 = c + z3

    e0, e1, e2, e3 = _tw_vecsum4(z00p, b0, b1, s3)
    r1, r2 = _tw_vseb2_3(e1, e2, e3)
    return e0, r1, r2
end

"""
    tw_prod33_fast_limbs(x::NTuple{3,Float64}, y::NTuple{3,Float64})

Fabiano–Muller–Picot Algorithm 10, fast triple-word × triple-word product.
The proven relative error bound (assuming the paper's TW input contract and no
under/overflow) is `44u^3 + 176u^4`.
"""
@inline function tw_prod33_fast_limbs(
    x::NTuple{3,Float64},
    y::NTuple{3,Float64},
)
    x0, x1, x2 = x
    y0, y1, y2 = y

    z00p, z00m = two_prod(x0, y0)
    z01p, z01m = two_prod(x0, y1)
    z10p, z10m = two_prod(x1, y0)

    b0, b1, b2 = _tw_vecsum3(z00m, z01p, z10p)
    c = fma(x1, y1, b2)

    z31 = fma(x0, y2, z10m)
    z32 = fma(x2, y0, z01m)
    z3 = z31 + z32
    s3 = c + z3

    e0, e1, e2, e3 = _tw_vecsum4(z00p, b0, b1, s3)
    r1, r2 = _tw_vseb2_3(e1, e2, e3)
    return e0, r1, r2
end

@inline function tw_prod23_fast(x::MultiFloat{Float64,2}, y::MultiFloat{Float64,3})
    return MultiFloat{Float64,3}(tw_prod23_fast_limbs(x._limbs, y._limbs))
end

@inline function tw_prod23_fast(x::MultiFloat{Float64,3}, y::MultiFloat{Float64,2})
    return tw_prod23_fast(y, x)
end

@inline function tw_prod33_fast(x::MultiFloat{Float64,3}, y::MultiFloat{Float64,3})
    return MultiFloat{Float64,3}(tw_prod33_fast_limbs(x._limbs, y._limbs))
end
