# Float64x5-x8 safe multiplication family

Status: **M4 Experimental correctness baseline**.

`Experimental.mul_safe` exists to establish a reviewable higher-limb
multiplication source of truth before fixed-cost or SIMD optimization. It is
intentionally allocating, sorting, and renormalizing far more components than a
future hot kernel should.

## Construction

For width `N in 5:8` and normalized finite Float64 inputs:

1. evaluate all `N^2` limb products using `MultiFloats.two_prod`;
2. retain both rounded product and residual from every pair (`2N^2` components);
3. canonicalize signed zeros to `+0.0`;
4. sort components by deterministic `(abs(value), value)` order, descending;
5. fully apply `MultiFloats.renormalize` to the exact component tuple;
6. retain the leading N normalized components;
7. renormalize that N-term head.

Width-specific wrappers `mul5_safe` ... `mul8_safe` call this same
implementation. No FastTwoSum appears in the family.

The deterministic sort makes `x*y` and `y*x` enter iterative renormalization with
the same component order, so bitwise commutativity is a permanent regression
gate rather than a tolerance-only expectation.

## Exactness gates

For every accepted test case, every TwoProd pair must satisfy exactly:

```text
Rational(p) + Rational(e) = Rational(x_i) * Rational(y_j)
```

The full product must then satisfy:

```text
value(returned N-limb head)
+ value(discarded normalized limbs N+1 : 2N^2)
= value(x) * value(y)
```

using exact `Rational{BigInt}` arithmetic.

Additional gates require:

- canonical component sequence identical under operand swap;
- full `2N^2` expansion normalized;
- returned N-limb expansion normalized;
- `mul_safe(x,y) === mul_safe(y,x)`;
- bitwise equality with independent M2 `reference_mul`;
- equality between generic and named width wrappers;
- exact zero/one/sign cases;
- dense multi-limb and scaled inputs;
- powers of two and near-boundary products;
- explicit domain errors outside N=5:8.

## Current domain

The baseline deliberately rejects pair computations with:

- nonfinite TwoProd product/residual components;
- nonzero subnormal rounded products;
- nonzero subnormal TwoProd residuals;
- nonzero input pairs whose rounded product underflows to zero.

These exclusions are conservative. Gradual-underflow behavior should be expanded
only after the pairwise exactness and final error contract are analyzed there.

## First width-specific error measurements

Zen 3 GitHub-hosted runner, Julia 1.10.11:

| Width | Dense/scaled cases | Max observed `|z-xy|/(u^N|xy|)` |
|---|---:|---:|
| x5 | 300 + 300 | 0.0484069 |
| x6 | 80 + 80 | 0.0183703 |
| x7 | 60 + 60 | 0.00725543 |
| x8 | 40 + 40 | 0.00253803 |

All measured dense/scaled cases had nonzero discarded tails, zero oracle
mismatches, zero normalization failures, and zero commutativity failures.
Power-of-two boundary corpora were exact.

CI requires each width's measured relative constant to remain below `C=1`. This
is a deliberately loose empirical regression gate, **not a worst-case theorem**.

## First timings

On the same runner:

- x5: 32.84 ms / 500 products;
- x6: 5.059 ms / 40 products;
- x7: 5.710 ms / 25 products;
- x8: 4.492 ms / 12 products.

The case counts differ deliberately. Per-operation cost grows quickly as the safe
expansion grows from 50 to 128 components. These measurements document a
correctness floor, not a performance target.

## Promotion rules

A future fixed-cost or SIMD multiplication candidate must preserve at minimum:

1. independent `reference_mul` agreement on permanent adversarial corpora;
2. normalized output and commutativity;
3. an explicit discarded-tail bound;
4. explicit overflow/subnormal/underflow semantics;
5. every FastTwoSum magnitude proof if FastTwoSum is introduced;
6. end-to-end downstream residual/certificate integrity before solver adoption.

The next arithmetic milestone is M5: correctness-first direct FMA x5, then x6-x8,
before any direct fused-network minimization.
