# M6 reciprocal, division, and square-root reference oracles

Status: **Experimental correctness oracle for Float32/Float64 x5-x8**.

These functions exist to validate future precision-doubling arithmetic. They are
not performance implementations.

## Why these oracles differ from M2

For dyadic MultiFloat inputs, add/mul/FMA have dyadic exact source results.
Reciprocal and division are generally non-dyadic, while square root is generally
irrational. Therefore a single finite BigFloat precision cannot honestly be
called an exact result.

M6 uses **precision-stabilized final packing** instead.

## API

Under `MultiFloatArithmetic.Experimental`:

```julia
reference_inv(x)
reference_div(x, y)
reference_sqrt(x)
```

The scalar functions support Float32/Float64 with `N in 5:8`. `MultiFloatVec`
wrappers are lane-wise calls to the scalar oracle so future vector Newton kernels
do not share implementation failure modes with their reference.

## Reciprocal and division

The input expansion is converted exactly to `Rational{BigInt}`.

- reciprocal source value: `1 / Rational(x)`;
- division source value: `Rational(x) / Rational(y)`.

The exact rational value is then evaluated in MPFR and packed into an N-limb
MultiFloat.

Because the rational is usually non-dyadic, packing starts at a conservative
width-dependent BigFloat precision and doubles precision until **two consecutive
N-limb packings are bitwise identical**. The implementation fails rather than
silently accepting a result that has not stabilized by 32768 bits.

## Square root

`reference_sqrt` converts the finite normalized input exactly to a rational and
evaluates

```text
sqrt(numerator / denominator)
```

in MPFR. It uses the same precision-doubling final-packing stabilization. Zero is
exact; negative inputs are rejected.

## Independent CI validation

The adaptive-stabilization result is not accepted solely because it agrees with
itself at two consecutive precisions. CI also computes a separate **16384-bit**
one-shot pack for reciprocal, division, and square root and requires bitwise
identity with the adaptive oracle.

The permanent test matrix includes:

- Float64x5, x6, x7, x8 scalar cases;
- Float32x5/x8 smoke cases;
- normalized-output checks;
- reciprocal/division/sqrt exact identities;
- Vec2 lane equivalence;
- zero-denominator and negative-sqrt domain errors;
- N<5/N>8 rejection;
- independence from ambient BigFloat precision.

The first M6 oracle implementation passes this suite on Linux Julia 1.10,
current Julia, and macOS current.

## Domain

Accepted scalar inputs are normalized finite Float32/Float64 expansions with
5-8 limbs.

- `reference_inv(0)` is rejected;
- `reference_div(x, 0)` is rejected;
- `reference_sqrt(x)` rejects negative x;
- sqrt(0) returns exact zero.

The reference layer does not yet make a performance or reproducibility claim for
candidate hardware division/sqrt instructions; it exists specifically to judge
those future candidates.

## Next M6 step

Develop **x8 reciprocal first** with a deliberate 1->2->4->8 refinement path:

1. Float64 seed;
2. x2/x4 refinement using trusted upstream arithmetic;
3. lift the approximation to x8;
4. final x8 Newton/Karp-Markstein corrections using accepted safe multiplication
   and direct-FMA/submul-style residuals;
5. compare against `reference_inv` at every accepted test case.

Only after reciprocal x8 is stable should division reuse reciprocal or a direct
quotient-correction path. Square root follows with its own residual formulation.
The candidate iteration must never be used as its own reference.
