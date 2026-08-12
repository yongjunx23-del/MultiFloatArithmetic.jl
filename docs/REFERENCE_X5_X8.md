# x5-x8 reference arithmetic contract

Status: **M2 correctness oracle**. This is not a performance API.

## Purpose

`MultiFloatArithmetic.Experimental.reference_add`, `reference_sub`,
`reference_mul`, and `reference_fma` provide an implementation-independent source
of truth for future 5- through 8-limb arithmetic research.

The reference path exists because an optimized candidate must be rejectable on
correctness before instruction count, SIMD throughput, or solver wall time is
considered. The x4 cancellation finding is the motivating example: broad random
tests and favorable benchmarks were not sufficient to validate a hidden
FastTwoSum ordering assumption.

## Domain

The current reference contract accepts:

- base limb type `Float32` or `Float64`;
- expansion length `N in 5:8`;
- finite inputs;
- inputs already normalized according to `MultiFloats.isnormalized`.

Other formats and widths are rejected explicitly rather than silently falling
through to an unvalidated path.

## Exact arithmetic model

Each scalar `MultiFloat` is converted with `Rational{BigInt}(x)`. Because every
binary floating-point limb is a dyadic rational, this conversion is exact.

The requested operation is then performed in rational arithmetic:

```text
add:  qx + qy
sub:  qx - qy
mul:  qx * qy
fma:  qx * qy + qc
```

Only after the exact operation is complete is the result converted back to an
N-limb `MultiFloat`.

## Packing precision

Packing uses an internally selected BigFloat precision

```text
p = 2 * span(T) + 2 * ndigits(N, base=2) + 32
```

where

```text
span(T) = exponent(floatmax(T)) - exponent(floatmin(T)) + precision(T).
```

`span(T)` covers the binary positions from the largest finite exponent through
the smallest subnormal bit. The factor of two covers product span; the remaining
terms provide carry/N margin.

This precision choice is itself tested rather than assumed. CI compares every
reference result in the x5-x8 scalar test matrix against an independently packed
8192-bit BigFloat result. The reference operation also sets its own BigFloat
precision, so callers cannot weaken it by changing ambient MPFR precision.

## Vector semantics

`MultiFloatVec` reference operations are defined lane by lane through the scalar
oracle. This is intentionally slow and avoids sharing a vector arithmetic
network with future SIMD candidates.

A vector candidate is accepted only if every lane agrees with the scalar oracle
under the candidate's stated error/rounding contract.

## Non-goals

The reference path does **not** aim to be:

- allocation free;
- branch free;
- SIMD optimized;
- suitable for production hot loops;
- a replacement for eventual verified `add5`-`add8`, `mul5`-`mul8`, or direct
  FMA networks.

Its cost is acceptable because it is a test/research oracle.

## Gate for future higher-limb kernels

Before a proposed x5-x8 arithmetic network can be promoted, it must at minimum:

1. match this oracle on exact identities and a permanent adversarial corpus;
2. produce normalized outputs;
3. satisfy an explicit error metric with a reviewable bound;
4. preserve required symmetries such as multiplication commutativity;
5. pass scalar/SIMD lane checks where applicable;
6. only then proceed to codegen and performance evaluation.
