# Float64x5-x8 safe direct FMA family

Status: **M5 Experimental correctness baseline**.

`Experimental.fma_safe` establishes direct higher-limb `x*y+c` behavior before
fixed-cost or SIMD fused optimization.

## Construction

For width `N in 5:8`:

1. retain all `2N^2` exact TwoProd product/residual components from the accepted
   safe multiplication decomposition;
2. append all N exact limbs of `c`;
3. canonicalize signed zero and deterministically sort the complete component
   multiset;
4. fully renormalize the `2N^2+N`-component direct expansion;
5. keep the leading N limbs and renormalize the head.

Component counts are 55, 78, 105, and 136 for x5-x8. No FastTwoSum is used.
The already-frozen `fma5_safe` wrapper remains as the x5 baseline; generic x5 is
required to match it bit-for-bit.

## Exactness gates

Every accepted case checks each TwoProd pair exactly, then requires

```text
value(returned N-limb head)
+ value(discarded normalized tail)
= value(x) * value(y) + value(c)
```

with exact `Rational{BigInt}` arithmetic.

Permanent gates also require normalized full/output expansions, bitwise
`reference_fma` equality, x/y symmetry, generic/wrapper equality, identities,
dense/scaled inputs, power-of-two cancellation, and destructive cancellation
past nominal N-limb precision.

## Current domain

The product side inherits the conservative `mul_safe` exclusions for nonfinite,
nonzero-subnormal, and underflow-to-zero TwoProd components. The addend must be
finite and normalized. These restrictions remain until gradual-underflow
behavior is explicitly proved.

## First empirical direct-FMA results

Zen 3, Julia 1.10.11:

| Width | Max observed `|err|/(u^N(|xy|+|c|))` | Direct oracle / norm / symmetry failures |
|---|---:|---:|
| x5 | 0.0492768 | 0 / 0 / 0 |
| x6 | 0.0227603 | 0 / 0 / 0 |
| x7 | 0.00452122 | 0 / 0 / 0 |
| x8 | 0.000825581 | 0 / 0 / 0 |

CI uses C=1 independently per width as an empirical regression gate. These are
not formal worst-case bounds.

## Why direct FMA is distinct from safe mul-then-add

The rounded composition

```julia
add_safe(mul_safe(x, y), c)
```

was compared with `reference_fma` on the same seeded diagnostics:

| Width | Ordinary | Scaled | Destructive cancellation |
|---|---:|---:|---:|
| x5 | 45 / 200 | 41 / 200 | 150 / 150 |
| x6 | 8 / 50 | 14 / 50 | 48 / 48 |
| x7 | 7 / 35 | 9 / 35 | 32 / 32 |
| x8 | 6 / 25 | 6 / 25 | 24 / 24 |

The composition disagreed with the exact direct-FMA oracle in **every tested
destructive-cancellation case at every width**, while direct `fma_safe` had zero
oracle mismatches. Separately correct rounded add/mul are therefore not an exact
substitute for direct FMA in residual, refinement, or certification paths.

## First timing A/B

Hosted-runner case counts differ by width because the exact expansion grows to
136 components at x8:

| Width | Direct safe FMA | Safe mul-then-add | composition/direct |
|---|---:|---:|---:|
| x5 | 9.852 ms / 120 | 7.993 ms / 120 | 0.811x |
| x6 | 3.044 ms / 20 | 2.550 ms / 20 | 0.837x |
| x7 | 3.186 ms / 12 | 2.726 ms / 12 | 0.856x |
| x8 | 2.625 ms / 6 | 2.314 ms / 6 | 0.882x |

The current direct correctness implementation is slower. Its present purpose is
preserving fused numerical behavior; fixed-cost optimization begins only after
formal direct-tail analysis.

## Promotion rules

A future direct-FMA performance candidate must preserve independent
`reference_fma` agreement, normalization, x/y symmetry, cancellation behavior,
explicit error/underflow semantics, and every FastTwoSum ordering proof if such
operations are introduced.
