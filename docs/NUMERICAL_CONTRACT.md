# Numerical contract

This document defines the acceptance boundary for the top-level
`MultiFloatArithmetic.fma_fast` kernel. It is intentionally narrower than the
set of values for which the Julia code happens to return a result.

## Arithmetic model

The research contract assumes:

- IEEE-style binary32 or binary64 arithmetic;
- round to nearest, ties to even;
- finite normalized `MultiFloat` input expansions;
- no compiler reassociation or `@fastmath` transformation;
- no flush-to-zero or denormals-are-zero mode in a path where gradual underflow
  matters.

NaN, infinity, signed-zero details, overflow, and underflow-adjacent behavior are
not currently part of the accepted contract. They require separate semantics
and tests before being promised.

## Input invariant

Each input expansion must satisfy the normalization/non-overlap invariant used
by `MultiFloats.jl`. The hot kernel does not revalidate this precondition because
runtime checks would change the branch-free network being measured.

## Error metric

For `z = fma_fast(x, y, c)`, the current empirical acceptance tests use an
operand-relative bound of the form

```text
|z - (x*y + c)| <= C_N * u^N * (|x*y| + |c|) + oracle_slack,
```

where `u` is the unit roundoff of the base format and the current test constants
are `C_2 = 34`, `C_3 = 184`, and `C_4 = 812`.

These constants are test gates, not yet machine-verified theorems. The project
must not describe them as formal worst-case guarantees until verifier artifacts
are checked into the repository.

## Destructive cancellation

When `x*y` and `c` nearly cancel, the exact result may be far smaller than
`|x*y| + |c|`. The operand-relative bound can therefore coexist with a large
result-relative error. This is expected behavior under the current contract.

Consequences:

- suitable target: throughput-oriented accumulation where the surrounding error
  analysis is operand-relative;
- unsuitable default: KKT residuals, iterative refinement corrections,
  stopping criteria, feasibility checks, and final certificates;
- required downstream policy: retain a stronger path for cancellation-sensitive
  calculations and compare both paths in original problem coordinates.

## SIMD semantics

`MultiFloatVec` evaluation is intended to be lane-wise equivalent to scalar
`fma_fast`. CI currently checks exact lane equality for widths 2, 4, and 8.
Other widths supported by upstream types are not yet part of the explicit test
matrix.

## Compiler and code-generation invariant

The source order of additions and error-free transforms is part of the
algorithm. Refactors that are algebraically equivalent over the reals are not
necessarily equivalent in floating-point arithmetic.

Every arithmetic-network change must therefore rerun:

1. deterministic and randomized correctness tests;
2. destructive-cancellation cases;
3. scalar/SIMD lane-equivalence tests;
4. native-code diagnostics;
5. architecture-specific performance measurements.

## Experimental operations

`Experimental.mul_scalar` and `Experimental.div_digits` are retained to
reproduce a completed A/B study. They are not performance recommendations. The
quotient-digit division candidate was slower than upstream division in every
Float64x2/x3/x4 scalar and Vec4 case measured on the 2026-08-11 hosted runner.

## Formal verification target

The next proof milestone is to encode the accumulation portion of each x2/x3/x4
network as a reviewable floating-point accumulation network and establish:

- normalized/non-overlapping outputs;
- explicit discarded-tail bounds;
- assumptions required by each `fast_two_sum` use;
- precision-independent validity for binary32/binary64 where claimed.

The repository should store verifier inputs, tool/version metadata, generated
proof summaries, and a CI command that reproduces them.
