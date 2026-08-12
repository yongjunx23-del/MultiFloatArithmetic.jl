# Numerical contract

This document defines the acceptance boundary for
`MultiFloatArithmetic.fma_fast`. It is intentionally narrower than the set of
values for which the Julia code happens to return a result.

## Arithmetic model

The research contract assumes:

- IEEE-style binary32 or binary64 arithmetic;
- round to nearest, ties to even;
- finite normalized `MultiFloat` input expansions;
- no compiler reassociation or `@fastmath` transformation;
- no flush-to-zero or denormals-are-zero mode where gradual underflow matters.

NaN, infinity, signed-zero details, overflow, and underflow-adjacent behavior are
not currently promised and need separate semantics/tests.

## Input invariant

Each input expansion must satisfy the normalization/non-overlap invariant used by
`MultiFloats.jl`. The arithmetic kernels do not revalidate this precondition in
the hot path.

## Error metric

For `z = fma_fast(x, y, c)`, empirical acceptance uses

```text
|z - (x*y + c)| <= C_N * u^N * (|x*y| + |c|) + oracle_slack,
```

with current test constants `C_2 = 34`, `C_3 = 184`, and `C_4 = 812`.

These are **test gates, not machine-proved worst-case theorems**. Formal
error-bound verification is still pending.

## Destructive cancellation

When `x*y` and `c` nearly cancel, the exact result can be much smaller than
`|x*y| + |c|`, so an operand-relative bound can coexist with a large
result-relative error. Consequently this API is not the default choice for KKT
residuals, iterative refinement, stopping criteria, feasibility checks, or final
certificates.

Cancellation also matters structurally. FPAN decomposition found that the
original x4 final compression used `fast_two_sum(b, a1)` without a universally
valid magnitude ordering. A seeded concrete Float64 corpus with `c ≈ -x*y`
showed 3,755/10,000 ordering violations and 5,258/10,000 non-normalized old
outputs.

## Per-width structural status

### x2

The fixed-cost x2 source network has a pinned FPAN structural proof for its
explicit FastTwoSum precondition and output non-overlap relation.

### x3

The original x3 network needed an additional final TwoSum/FastTwoSum pair. The
repaired fixed-cost x3 network has pinned FPAN proofs for the final non-overlap
relations.

### x4

The accepted correctness baseline does **not** rely on the invalid old
FastTwoSum ordering. Its final compression uses general TwoSum transforms and
then calls `MultiFloats.renormalize`.

This means x4 is currently different from x2/x3:

- correctness takes priority over fixed-cost branch-free execution;
- SIMD remains useful empirically, but scalar x4 is slower than upstream;
- a future fixed-cost x4 replacement must prove source-specific intermediate
  bounds or use an established expansion-compression theorem;
- assumption-free short cleanup networks tested at p=24 and p=53 did not prove
  the top two output non-overlap relations.

## SIMD semantics

`MultiFloatVec` evaluation is intended to be lane-wise equivalent to scalar
`fma_fast`. CI checks exact lane equality for widths 2, 4, and 8. The x4
`renormalize` fallback uses upstream vector renormalization semantics.

## Compiler/code-generation invariant

Source ordering of additions and error-free transforms is part of the algorithm.
Algebraically equivalent refactors can be numerically different.

Every arithmetic-network change must rerun:

1. deterministic and randomized correctness tests;
2. destructive-cancellation cases;
3. scalar/SIMD lane-equivalence tests;
4. code-generation diagnostics;
5. architecture-specific performance measurements.

## Experimental operations

`Experimental.mul_scalar` and `Experimental.div_digits` are retained only for
reproducibility. Quotient-digit scalar division is not competitive and Vec4
results are architecture/type dependent.

## Remaining formal target

The next proof milestone is to establish reviewable discarded-tail/error bounds
for x2/x3, and source-specific intermediate bounds for a fixed-cost x4
compression. Verifier inputs, tool/version metadata, and proof summaries should
remain checked into the repository.
