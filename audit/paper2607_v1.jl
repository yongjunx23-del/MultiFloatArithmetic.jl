module Paper2607V1Audit

using MultiFloats
import MultiFloats: MultiFloat, MultiFloatVec, fast_two_sum, two_prod, two_sum

export paper_dw_fma, paper_dw_fma_limbs
export paper_qw_fma, paper_qw_fma_limbs
export paper_qw_fast_preconditions

# arXiv:2607.11391v1 Algorithm 1 / Appendix B.1: DW FMA, 17 flops.
@inline function paper_dw_fma_limbs(
    x::NTuple{2,T}, y::NTuple{2,T}, c::NTuple{2,T},
) where {T}
    x0, x1 = x
    y0, y1 = y
    c0, c1 = c

    p00, e00 = two_prod(x0, y0)
    p01 = x0 * y1
    p10 = x1 * y0
    ell = p01 + p10
    v = e00 + c1
    w = v + ell
    s, t = two_sum(p00, c0)
    tp = t + w
    return fast_two_sum(s, tp)
end

@inline paper_dw_fma(x::MultiFloat{T,2}, y::MultiFloat{T,2}, c::MultiFloat{T,2}) where {T} =
    MultiFloat{T,2}(paper_dw_fma_limbs(x._limbs, y._limbs, c._limbs))

@inline paper_dw_fma(
    x::MultiFloatVec{W,T,2}, y::MultiFloatVec{W,T,2}, c::MultiFloatVec{W,T,2},
) where {W,T} = MultiFloatVec{W,T,2}(paper_dw_fma_limbs(x._limbs, y._limbs, c._limbs))

# arXiv:2607.11391v1 Algorithm 3 / Appendix B.1: symmetric QW FMA, 146 flops.
# This is intentionally reproduced as published, including the single QW
# renormalization chain. The paper itself notes in Sec. 2.9.5 that the leading
# output pair need not be non-overlapping under extreme cancellation.
@inline function paper_qw_fma_limbs(
    x::NTuple{4,T}, y::NTuple{4,T}, c::NTuple{4,T},
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
    d = (p03 + p30) + (p12 + p21)

    b, r = two_sum(p00, c0)

    a1, f1 = two_sum(p01, p10)
    a1, f2 = two_sum(a1, e00)
    a1, f3 = two_sum(a1, c1)
    a1, f4 = two_sum(a1, r)

    a2, g1 = two_sum(p02, p20)
    a2, g2 = two_sum(a2, p11)
    etilde, g4 = two_sum(e01, e10)
    a2, g3 = two_sum(a2, etilde)
    a2, g5 = two_sum(a2, c2)
    a2, g6 = two_sum(a2, f1)
    a2, g7 = two_sum(a2, f2)
    a2, g8 = two_sum(a2, f3)
    a2, g9 = two_sum(a2, f4)

    t1 = e02 + e20
    t2 = e11 + d
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

    w0, w1 = fast_two_sum(b, a1)
    w1, w2 = two_sum(w1, a2)
    w2, w3 = two_sum(w2, a3)
    z0, rho = two_sum(w0, w1)
    z1, sig = two_sum(rho, w2)
    z2, z3 = fast_two_sum(sig, w3)
    return z0, z1, z2, z3
end

@inline paper_qw_fma(x::MultiFloat{T,4}, y::MultiFloat{T,4}, c::MultiFloat{T,4}) where {T} =
    MultiFloat{T,4}(paper_qw_fma_limbs(x._limbs, y._limbs, c._limbs))

@inline paper_qw_fma(
    x::MultiFloatVec{W,T,4}, y::MultiFloatVec{W,T,4}, c::MultiFloatVec{W,T,4},
) where {W,T} = MultiFloatVec{W,T,4}(paper_qw_fma_limbs(x._limbs, y._limbs, c._limbs))

# The paper/FPANVerifier FastTwoSum condition is exponent(a) >= exponent(b),
# or either operand equal to zero. It is not the stronger textbook sufficient
# check |a| >= |b| that an earlier project diagnostic used.
@inline function _paper_f2s_ok(a::Float64, b::Float64)
    (iszero(a) || iszero(b)) && return true
    return exponent(a) >= exponent(b)
end

# Replay exactly the two QW FastTwoSum gates from Algorithm 3 and return their
# concrete exponent-contract results. This is a falsification diagnostic only;
# the universal proof is performed independently by FPANVerifier.
function paper_qw_fast_preconditions(
    x::MultiFloat{Float64,4}, y::MultiFloat{Float64,4}, c::MultiFloat{Float64,4},
)
    x0, x1, x2, x3 = x._limbs
    y0, y1, y2, y3 = y._limbs
    c0, c1, c2, c3 = c._limbs

    p00, e00 = two_prod(x0, y0)
    p01, e01 = two_prod(x0, y1)
    p10, e10 = two_prod(x1, y0)
    p02, e02 = two_prod(x0, y2)
    p11, e11 = two_prod(x1, y1)
    p20, e20 = two_prod(x2, y0)
    d = ((x0*y3) + (x3*y0)) + ((x1*y2) + (x2*y1))

    b, r = two_sum(p00, c0)
    a1, f1 = two_sum(p01, p10)
    a1, f2 = two_sum(a1, e00)
    a1, f3 = two_sum(a1, c1)
    a1, f4 = two_sum(a1, r)

    a2, g1 = two_sum(p02, p20)
    a2, g2 = two_sum(a2, p11)
    etilde, g4 = two_sum(e01, e10)
    a2, g3 = two_sum(a2, etilde)
    a2, g5 = two_sum(a2, c2)
    a2, g6 = two_sum(a2, f1)
    a2, g7 = two_sum(a2, f2)
    a2, g8 = two_sum(a2, f3)
    a2, g9 = two_sum(a2, f4)

    t1 = e02 + e20
    t2 = e11 + d
    t3 = (t1 + t2) + c3
    t1 = g1 + g2
    t2 = g3 + g4
    t1 += t2
    t2 = g6 + g7
    t4 = g8 + g9
    t2 += t4
    t1 += t2
    t1 += g5
    a3 = t3 + t1

    first_ok = _paper_f2s_ok(b, a1)
    w0, w1 = fast_two_sum(b, a1)
    w1, w2 = two_sum(w1, a2)
    w2, w3 = two_sum(w2, a3)
    _, rho = two_sum(w0, w1)
    _, sig = two_sum(rho, w2)
    last_ok = _paper_f2s_ok(sig, w3)
    return first_ok, last_ok
end

end # module
