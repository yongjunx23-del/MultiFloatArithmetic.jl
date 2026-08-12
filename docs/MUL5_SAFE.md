# Float64x5 safe multiplication baseline

Status: **M4 Experimental correctness baseline**.

`Experimental.mul5_safe` exists to establish a defensible higher-limb
multiplication source of truth before fixed-cost or SIMD optimization. It is
intentionally slow and branchy.

## Construction

For normalized finite Float64x5 inputs `x` and `y`:

1. evaluate all 25 `x[i] * y[j]` limb pairs using `MultiFloats.two_prod`;
2. retain both the rounded product and residual from every pair, producing 50
   components;
3. canonicalize signed zero to `+0.0`;
4. sort the 50 finite components by deterministic `(abs(value), value)` order,
   descending;
5. fully apply `MultiFloats.renormalize` to the 50-term expansion;
6. retain the leading five normalized components;
7. renormalize that five-term head once more.

The sort is not a performance design. It gives `x*y` and `y*x` the same
component sequence before iterative renormalization, so commutativity can be a
bitwise regression gate rather than a tolerance-only property.

No FastTwoSum appears in this construction.

## Exactness gates

The test suite checks each accepted pair product individually:

```text
Rational(p) + Rational(e) == Rational(x_i) * Rational(y_j)
```

for every TwoProd `(p,e)` in the 5x5 grid.

It then checks the full multiplication identity:

```text
value(returned five-limb head)
+ value(discarded normalized limbs 6:50)
= value(x) * value(y)
```

with exact `Rational{BigInt}` arithmetic.

Additional permanent gates require:

- full 50-term expansion normalized;
- returned five-limb expansion normalized;
- canonical 50-term sequence identical under operand swap;
- `mul5_safe(x,y) === mul5_safe(y,x)`;
- bitwise equality with independent M2 `reference_mul`;
- exact zero and one identities;
- dense multi-limb inputs;
- exponent-scaled inputs;
- sign changes;
- powers of two and nearby five-limb boundaries;
- explicit domain errors for invalid inputs.

## Current domain

The baseline accepts normalized finite Float64x5 inputs only.

It deliberately rejects pair computations with:

- nonfinite TwoProd product/residual components;
- nonzero subnormal rounded products;
- nonzero subnormal TwoProd residuals;
- nonzero input pairs whose rounded product underflows to zero.

These exclusions are conservative. They avoid silently treating gradual
underflow as covered before the exact TwoProd decomposition and final-error
behavior are formally analyzed there.

## First empirical error measurement

Zen 3 GitHub-hosted runner, Julia 1.10.11:

| Corpus | Cases | Max `|z-xy|/(u^5|xy|)` |
|---|---:|---:|
| dense ordinary | 300 | 0.0484069 |
| scaled | 300 | 0.0417276 |
| powers of two | 25 | 0 |

All three corpora had:

- zero oracle mismatches;
- zero normalization failures;
- zero commutativity failures.

Every dense/scaled case had a nonzero discarded tail, so the oracle agreement is
not an artifact of exact representability.

CI currently requires the measured relative constant to remain below `C=1`.
That is a deliberately loose empirical regression gate, **not a worst-case
proof**.

## First timing

The same runner measured approximately:

```text
33.557 ms / 500 Float64x5 products
```

This is not competitive and is not meant to be. The current algorithm allocates,
sorts 50 terms, and repeatedly renormalizes them. Its purpose is to give M4 a
reviewable numeric baseline from which product components can later be removed or
reordered one change at a time.

## Promotion rules

Do not promote or optimize a replacement merely because it is faster. A future
Float64x5 multiplication candidate must at minimum preserve:

1. exact-oracle agreement on the permanent corpus;
2. normalized output;
3. commutativity;
4. a reviewable discarded-tail bound;
5. explicit underflow/overflow semantics;
6. every FastTwoSum ordering proof if FastTwoSum is introduced.

The next M4 step is safe x6-x8 multiplication with width-specific diagnostics,
then a formal multiplication tail bound before fixed-cost search.
