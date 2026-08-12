using MultiFloatArithmetic
using MultiFloats
using Random

const E = MultiFloatArithmetic.Experimental
const T = MultiFloat{Float64,5}

function pow2d(e::Int)
    e >= 0 ? (BigInt(1) << e) // BigInt(1) : BigInt(1) // (BigInt(1) << (-e))
end

function dense(; words=5, emin=-80, emax=80)
    numerator = BigInt(0)
    for _ in 1:words
        numerator = (numerator << 64) + BigInt(rand(UInt64))
    end
    numerator |= BigInt(1) << (64words - 1)
    q = numerator // (BigInt(1) << 64words)
    q *= pow2d(rand(emin:emax))
    q = rand(Bool) ? q : -q
    return setprecision(BigFloat, 1536) do
        T(BigFloat(q))
    end
end

function main()
    Random.seed!(0x6d15_2026)
    one_fail = 0
    self_fail = 0
    zero_fail = 0
    negleft_fail = 0
    negright_fail = 0
    ref_identity_fail = 0

    o = one(T)
    z = zero(T)
    for _ in 1:30
        x = dense()
        nx = setprecision(BigFloat, 1536) do
            T(BigFloat(-Rational{BigInt}(x)))
        end
        q_one = E.div5_safe(x, o)
        q_self = E.div5_safe(x, x)
        q_zero = E.div5_safe(z, x)
        q_nl = E.div5_safe(nx, x)
        q_nr = E.div5_safe(x, nx)
        one_fail += q_one !== x
        self_fail += q_self !== o
        zero_fail += q_zero !== z
        negleft_fail += q_nl != -o
        negright_fail += q_nr != -o
        ref_identity_fail += q_one !== E.reference_div(x, o)
        ref_identity_fail += q_self !== E.reference_div(x, x)
        ref_identity_fail += q_zero !== E.reference_div(z, x)
        ref_identity_fail += q_nl !== E.reference_div(nx, x)
        ref_identity_fail += q_nr !== E.reference_div(x, nx)
    end

    power_fail = 0
    power_ref_fail = 0
    for ex in (-120, -60, 0, 60, 120), ey in (-120, -60, 0, 60, 120)
        x = T(BigFloat(pow2d(ex)))
        y = T(BigFloat(pow2d(ey)))
        q = E.div5_safe(x, y)
        expected = T(BigFloat(pow2d(ex - ey)))
        power_fail += q !== expected
        power_ref_fail += q !== E.reference_div(x, y)
    end

    println("M6 Float64x5 division edge diagnostic")
    println("x/1 strict failures: ", one_fail)
    println("x/x strict failures: ", self_fail)
    println("0/x strict failures: ", zero_fail)
    println("(-x)/x numerical failures: ", negleft_fail)
    println("x/(-x) numerical failures: ", negright_fail)
    println("identity cases differing from reference_div: ", ref_identity_fail)
    println("power-of-two strict expected failures: ", power_fail)
    println("power-of-two cases differing from reference_div: ", power_ref_fail)

    # Domain classification is printed explicitly, not asserted here.
    cases = [
        ("zero denominator", () -> E.div5_safe(T(1.0), T(0.0))),
        ("infinite numerator", () -> E.div5_safe(T(Inf), T(1.0))),
        ("infinite denominator", () -> E.div5_safe(T(1.0), T(Inf))),
        ("subnormal denominator", () -> E.div5_safe(T(1.0), T(nextfloat(0.0)))),
        ("floatmax denominator", () -> E.div5_safe(T(1.0), T(floatmax(Float64)))),
    ]
    for (label, f) in cases
        outcome = try
            f()
            "returned"
        catch err
            string(typeof(err))
        end
        println(label, ": ", outcome)
    end
end

main()
