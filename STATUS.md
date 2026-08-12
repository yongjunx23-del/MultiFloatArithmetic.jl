# Project status

Status date: **2026-08-12**
Current milestone: **M3 — verified-safe higher-limb addition baselines**

## Executive assessment

M1.1 and M2 are merged. The project now has:

- pinned structural verification for the fixed-cost x2/x3 FMA networks;
- a cancellation-safe x4 FMA baseline after the original FastTwoSum defect was
  reproduced on concrete Float64 inputs;
- exact-rational x5-x8 Experimental add/sub/mul/FMA oracles;
- a correctness-first Float64x5 addition baseline whose returned value matches
  the exact M2 five-limb oracle on the permanent corpus.

Optimized x5-x8 arithmetic, formal error proofs, and native linear algebra are
still future milestones.

## M3 Float64x5 safe addition

`Experimental.add5_safe` deliberately avoids FastTwoSum. It:

1. applies `two_sum` to each of the five same-order limb pairs;
2. fully renormalizes the resulting exact ten-term expansion;
3. retains the leading five terms;
4. renormalizes the five-term head.

The low five normalized terms remain available to diagnostics, so CI checks the
exact dyadic identity

```text
value(result) + value(discarded_tail) = value(x) + value(y).
```

The permanent corpus includes signed random inputs, wide exponent separation,
exact identities, deep cancellation, power-of-two/carry boundaries, and highly
unbalanced operands. It requires normalized output, bitwise commutativity, and
bitwise equality with `Experimental.reference_add`.

### First discarded-tail measurement

On the initial 2,500-case diagnostic corpus:

- oracle mismatches: **0**;
- normalization failures: **0**;
- commutativity failures: **0**;
- ordinary cases with nonzero discarded tail: **583 / 1,000**;
- wide-exponent cases with nonzero discarded tail: **968 / 1,000**;
- near-cancellation cases with nonzero discarded tail: **0 / 500**;
- maximum observed `|err| / (u^5 (|x| + |y|))`: **~0.05284**.

CI uses `C = 1` only as a conservative empirical regression gate. It is not yet
a theorem. See `docs/ADD5_SAFE.md`.

The safe implementation measured roughly 0.5 ms for 5,000 scalar additions on
the first Zen 3 hosted-runner snapshot. Performance is not yet an acceptance
criterion for this baseline.

## M2 exact oracle

`MultiFloatArithmetic.Experimental` provides `reference_add`, `reference_sub`,
`reference_mul`, and `reference_fma` for Float32/Float64 expansions with
`N in 5:8`.

Each scalar input is converted to exact `Rational{BigInt}` dyadics, the operation
is performed exactly, and the result is packed once. CI cross-checks the packing
against an independent 8192-bit BigFloat result. Vector reference methods are
lane-wise scalar calls to avoid sharing implementation failure modes with future
SIMD candidates.

See `docs/REFERENCE_X5_X8.md`.

## x4 finding retained as a design constraint

The original x4 final compression began with `fast_two_sum(b, a1)`. FPAN
decomposition refuted the required ordering under cancellation. In 10,000 seeded
cases with `c` chosen as the x4 representation of `-x*y`:

- first-FastTwoSum ordering violations: **3,755 / 10,000**;
- old x4 non-normalized outputs: **5,258 / 10,000**;
- safe-baseline non-normalized outputs: **0 / 10,000**.

Higher-limb candidates may not inherit FastTwoSum assumptions without proving
the source-specific magnitude ordering.

## Current decisions

| Component | Decision |
|---|---|
| x2 scalar `fma_fast` | marginal architecture-dependent opt-in candidate |
| x3 scalar `fma_fast` | structurally verified but slower; do not auto-select |
| x2/x3 SIMD `fma_fast` | continue architecture-specific optimization path |
| x4 scalar `fma_fast` | correctness baseline only; upstream is faster |
| x4 SIMD `fma_fast` | safe baseline retains useful speedup |
| x5-x8 `reference_*` | exact Experimental correctness oracle |
| `add5_safe` | accepted M3 Experimental correctness baseline; not a fast kernel |
| add6-add8 safe baselines | next M3 step |
| optimized x5-x8 multiplication/FMA | not started |
| quotient-digit division | rejected as default; retain under `Experimental` |
| native linear algebra backend | not implemented |

## Acceptance gates for higher-limb kernels

1. domain, rounding, underflow, normalization, and error metric are explicit;
2. permanent adversarial cases agree with the exact M2 oracle;
3. normalized output and required symmetries pass;
4. every FastTwoSum precondition is proved or FastTwoSum is avoided;
5. discarded-tail/error constants are measured and then formally bounded;
6. SIMD lane semantics pass where applicable;
7. only then inspect codegen and performance;
8. downstream solver A/B must not degrade residuals or certificates.

## Immediate next milestone

Freeze the safe x5 baseline, then generalize the same no-FastTwoSum construction
to x6, x7, and x8. Each width must independently match `reference_add`, preserve
normalization/commutativity under the adversarial corpus, and establish an
empirical discarded-tail gate before any fixed-cost minimization begins.
