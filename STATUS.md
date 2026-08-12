# Project status

Status date: **2026-08-12**
Current milestone: **M5 — safe Float64x5 direct FMA baseline accepted**

## Executive assessment

The project now has correctness-first higher-limb baselines for addition and
multiplication through Float64x8, and the first direct higher-limb FMA baseline at
Float64x5.

- x2/x3 fixed-cost FMA networks have pinned structural verification;
- x4 uses a cancellation-safe TwoSum + `renormalize` FMA baseline;
- x5-x8 have exact-rational Experimental add/sub/mul/FMA oracles;
- x5-x8 have common no-FastTwoSum safe addition and multiplication families;
- Float64x5 now has a direct 55-component safe FMA baseline that combines exact
  product components with the addend before any five-limb truncation.

Formal high-limb tail proofs, fixed-cost x5-x8 add/mul/FMA, and native linear
algebra remain future work.

## M5 Float64x5 direct FMA

`Experimental.fma5_safe` does **not** implement FMA as `mul_safe` followed by
`add_safe`. Instead it:

1. obtains the exact 50 product/residual components from the 25 x5 TwoProd pairs;
2. appends all five exact limbs of `c`;
3. canonicalizes signed zeros and sorts the 55 components deterministically;
4. fully renormalizes the combined exact `x*y+c` expansion;
5. retains the leading five limbs and renormalizes the head.

No FastTwoSum is used. Product-side overflow/subnormal/underflow exclusions
inherit the current M4 `mul_safe` contract.

The permanent corpus checks each TwoProd pair exactly and requires:

```text
value(result) + value(discarded 50 normalized limbs)
= value(x) * value(y) + value(c)
```

It also requires normalized full/output expansions, bitwise equality with
`reference_fma`, bitwise x/y symmetry, exact identities, dense/scaled inputs,
power boundaries, and destructive cancellation beyond the nominal five-limb
scale.

### First direct-FMA diagnostic

Zen 3, Julia 1.10.11:

| Corpus | Cases | Direct oracle failures | Mul-then-add oracle failures | Max `|err|/(u^5(|xy|+|c|))` |
|---|---:|---:|---:|---:|
| dense ordinary | 200 | 0 | 45 | 0.0456757 |
| scaled | 200 | 0 | 41 | 0.0492768 |
| destructive cancellation | 150 | 0 | 150 | ~1.0e-67 |

The direct path had zero normalization and x/y-symmetry failures. Its maximum
observed result-relative error was of order `1e-81` in these corpora.

The composition result is the important M5 design signal: even though both
`mul_safe` and `add_safe` independently match their operation-specific exact
oracles, the intermediate five-limb multiplication rounding causes the composed
`add_safe(mul_safe(x,y),c)` path to differ from the exact direct-FMA oracle in
**100% of the destructive-cancellation diagnostic**.

CI uses `C=1` for the direct operand-relative metric as an empirical regression
gate, not a theorem.

### First timing A/B

For 120 Float64x5 cases on the same runner:

- direct 55-component safe FMA: **7.707 ms**;
- safe multiplication followed by safe addition: **6.231 ms**;
- composed/direct = **0.809x**.

Thus the current direct correctness baseline is about 24% slower than the
already-slow safe composition. M5 therefore prioritizes removal of the
intermediate-rounding error first; fixed-cost direct-FMA optimization comes only
after the safe family and formal tail analysis exist.

## M4 safe multiplication family

`Experimental.mul_safe` and `mul5_safe` ... `mul8_safe` keep all `2N^2` TwoProd
components, canonicalize their order, fully renormalize, and truncate only after
the exact product expansion is formed. No FastTwoSum is used.

First max observed `|err|/(u^N|xy|)` values were 0.0484069 (x5), 0.0183703
(x6), 0.00725543 (x7), and 0.00253803 (x8), with zero measured oracle,
normalization, or commutativity failures.

## M3 safe addition family

`Experimental.add_safe` and `add5_safe` ... `add8_safe` use N general TwoSums,
full 2N-term exact renormalization, N-limb truncation, and head renormalization.
First max observed `|err|/(u^N(|x|+|y|))` values were 0.0528336, 0.0254805,
0.0109878, and 0.00541727 for x5-x8.

## Current decisions

| Component | Decision |
|---|---|
| x2/x3 fixed-cost FMA | structurally verified; architecture-sensitive performance |
| x4 FMA | cancellation-safe correctness baseline |
| x5-x8 `reference_*` | exact Experimental oracle |
| x5-x8 `add_safe` | accepted M3 Experimental correctness family |
| x5-x8 `mul_safe` | accepted M4 Experimental correctness family |
| `fma5_safe` | accepted M5 Experimental direct-FMA correctness baseline |
| mul-safe -> add-safe FMA composition | not an exact direct-FMA substitute; fails all tested destructive-cancellation oracle cases |
| fixed-cost high-limb add/mul/FMA | formal tail analysis first |
| safe direct FMA x6-x8 | next M5 step |
| reciprocal/div/sqrt | M6 |
| native DOT/SYRK/TRSM/GEMM | M7 |

## Immediate next milestone

Generalize the direct-FMA construction to Float64x6, x7, and x8. Width N should
combine all `2N^2` exact TwoProd product/residual components with N addend limbs,
forming `2N^2+N` components (78, 105, and 136 for x6-x8), then fully normalize
before N-limb truncation. Each width must independently pass exact head+tail,
`reference_fma`, x/y-symmetry, destructive-cancellation, and empirical
`u^N(|xy|+|c|)` gates before any fused-network minimization.
