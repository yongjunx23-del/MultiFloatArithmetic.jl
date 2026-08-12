# Float64x5 safe direct FMA baseline

Status: **M5 Experimental correctness baseline**.

`Experimental.fma5_safe` establishes a direct higher-limb `x*y+c` source of
truth before fixed-cost or SIMD fused optimization. It is intentionally branchy,
allocating, and slower than the separately rounded safe composition.

## Direct construction

For normalized finite Float64x5 inputs `x`, `y`, and `c`:

1. evaluate all 25 `x[i]*y[j]` pairs with the M4 exact TwoProd decomposition;
2. retain all 25 rounded products and all 25 residuals (50 product components);
3. append all five exact limbs of `c`;
4. canonicalize signed zero and sort all 55 components by deterministic
   `(abs(value), value)` order;
5. fully `renormalize` the combined exact `x*y+c` expansion;
6. retain the leading five normalized components;
7. renormalize that five-term head.

The algorithm reuses only `_mul_safe_terms(x,y)`, not the rounded five-limb
`mul_safe(x,y)` result. Thus multiplication is not truncated before the addend is
introduced.

No FastTwoSum appears in the construction.

## Exactness gates

Every accepted case first checks every TwoProd pair exactly:

```text
Rational(p) + Rational(e) = Rational(x_i) * Rational(y_j)
```

Then the direct FMA must satisfy:

```text
value(returned five-limb head)
+ value(discarded normalized limbs 6:55)
= value(x) * value(y) + value(c)
```

with exact `Rational{BigInt}` arithmetic.

Permanent gates additionally require:

- 55-component sequence identical under `x <-> y`;
- full 55-term expansion normalized;
- returned five-limb result normalized;
- bitwise equality with independent `reference_fma`;
- bitwise x/y symmetry;
- zero/one/addend identities;
- dense and exponent-scaled inputs;
- power-of-two/boundary cases;
- destructive cancellation to depths beyond the nominal five-limb scale;
- explicit invalid/overflow/underflow-domain errors.

## Current domain

The product component decomposition inherits M4's conservative exclusions:
nonfinite TwoProd components, nonzero subnormal products/residuals, and pair
products that underflow to zero are rejected pending dedicated gradual-underflow
analysis. The addend must be finite and normalized.

## First empirical error measurement

Zen 3, Julia 1.10.11:

| Corpus | Cases | Direct oracle failures | Max `|z-(xy+c)|/(u^5(|xy|+|c|))` |
|---|---:|---:|---:|
| dense ordinary | 200 | 0 | 0.0456757 |
| scaled | 200 | 0 | 0.0492768 |
| destructive cancellation | 150 | 0 | ~1.01e-67 |

There were zero normalization and x/y-symmetry failures. Direct result-relative
errors in the seeded corpora were of order `1e-81`.

CI currently requires the operand-relative constant to stay below `C=1`. This is
a deliberately loose empirical regression gate, **not a worst-case theorem**.

## Direct versus rounded mul-then-add

The same diagnostics also computed:

```julia
add_safe(mul_safe(x, y), c)
```

and compared it with `reference_fma`.

| Corpus | Composition oracle mismatches |
|---|---:|
| ordinary | 45 / 200 |
| scaled | 41 / 200 |
| destructive cancellation | 150 / 150 |

This is the key numerical justification for M5. `mul_safe` and `add_safe` are
individually oracle-matching operation baselines, but the five-limb rounding
between them changes the FMA result. Under the seeded destructive-cancellation
corpus it changed **every** result relative to the direct exact-FMA oracle.

For residual/refinement/certification paths, a direct FMA therefore cannot be
replaced by the rounded composition merely because both component operations are
correct in isolation.

## First timing

For 120 cases on the same runner:

```text
direct 55-component FMA: 7.707 ms
safe mul-then-add:        6.231 ms
composition/direct:       0.809x
```

The direct correctness baseline is roughly 24% slower. This is expected and not
a rejection: the next optimization problem is to remove unnecessary product and
normalization work while preserving the direct oracle behavior.

## Promotion rules

A future x5 direct-FMA performance candidate must preserve:

1. independent `reference_fma` agreement on permanent adversarial/cancellation
   corpora;
2. normalized output and x/y symmetry;
3. explicit discarded-tail and cancellation-aware error bounds;
4. explicit overflow/subnormal/underflow semantics;
5. every FastTwoSum ordering proof if FastTwoSum is introduced;
6. downstream residual/certificate integrity before solver adoption.

The next M5 step is the same safe direct construction for x6-x8, then a formal
direct-FMA tail analysis before fixed-cost fused-network search.
