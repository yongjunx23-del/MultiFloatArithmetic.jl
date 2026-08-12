# Verified Float64x5-Float64x8 arithmetic research plan

Research-only until formal verification and numerical gates pass.

## Core rule

The x4 cancellation defect permanently sets the search order: correctness,
normalization, exact-oracle agreement, and explicit error accounting before
branch-free/fixed-cost minimization. Every FastTwoSum must have a proved
source-specific magnitude precondition.

## M2 exact oracle — complete baseline

`Experimental.reference_add/sub/mul/fma` provides the independent x5-x8 rejection
standard via exact `Rational{BigInt}` arithmetic and one final pack.

## M3 addition — safe family complete baseline

`Experimental.add_safe` covers Float64 x5-x8 with exact 2N-term accumulation and
full renormalization before N-limb truncation. First empirical constants in
`u^N(|x|+|y|)` units are ~0.05283, 0.02548, 0.01099, and 0.00542 for x5-x8.

## M4 multiplication — safe family complete baseline

`Experimental.mul_safe` covers Float64 x5-x8 by preserving all `2N^2` TwoProd
components, canonicalizing their order, fully renormalizing, then truncating.
First empirical relative constants are ~0.04841, 0.01837, 0.00726, and 0.00254.

The x8 construction preserves/sorts/renormalizes 128 components. It is a
correctness reference, not a performance design: the first requested comparison
measured it roughly 3,600x slower than BigFloat512 multiplication on the hosted
Zen 3 runner.

## M5 direct FMA — safe x5-x8 family complete baseline

`Experimental.fma_safe` deliberately does not call the rounded N-limb
multiplication result. For width N it combines:

- all `2N^2` exact product/residual components from the M4 TwoProd decomposition;
- all N exact limbs of `c`;
- deterministic canonical ordering;
- one full renormalization of the resulting direct expansion;
- N-limb truncation only after the exact `x*y+c` expansion exists.

Full component counts are 55, 78, 105, and 136 for x5-x8.

Permanent gates require exact TwoProd reconstruction, exact direct-FMA head+tail
accounting, normalized full/output expansions, `reference_fma` equality, x/y
symmetry, generic/wrapper equality, identities, dense/scaled cases, powers of
two, and deep destructive cancellation.

First maximum direct-FMA `|err|/(u^N(|xy|+|c|))` values over the ordinary/scaled
seeded corpora:

- x5: ~0.0492768
- x6: ~0.0227603
- x7: ~0.00452122
- x8: ~0.000825581

All widths had zero direct oracle, normalization, and x/y-symmetry failures.

### Why direct FMA is mandatory for later residual/refinement work

The rounded composition `add_safe(mul_safe(x,y),c)` disagreed with the exact
`reference_fma` oracle in every destructive-cancellation case tested:

- x5: 150 / 150
- x6: 48 / 48
- x7: 32 / 32
- x8: 24 / 24

The direct path had zero oracle mismatches. Therefore separately correct add/mul
baselines are not a substitute for direct FMA in cancellation-sensitive residual,
refinement, stopping-criterion, or certificate code.

The safe direct path remains slower than safe composition. First
composition/direct timing ratios were 0.811x, 0.837x, 0.856x, and 0.882x for
x5-x8. This is acceptable only because M5 is a rejection/reference baseline.

## Post-M5 fixed-cost research track

Now that the safe families exist, fixed-cost search can proceed without losing a
correctness reference. The priority order is:

1. derive reviewable tail/error bounds for the safe add/mul/direct-FMA families;
2. reduce x8 multiplication first because it is the dominant measured
   performance failure versus BigFloat512;
3. search diagonal/staged product networks that retain exact-accounting hooks;
4. prove every proposed FastTwoSum ordering from source-specific bounds;
5. compare each minimized network bitwise against the safe/reference oracles on
   the permanent adversarial corpus;
6. only then introduce SIMD/MultiFloatVec variants and codegen gates.

A fixed-cost candidate is rejected immediately if it loses normalization,
oracle equality, required symmetry, cancellation behavior, or its stated error
bound, regardless of instruction count.

## M6 — reciprocal/division/sqrt correctness baseline next

Start with Float64x5 reciprocal and division before generalizing widths.
Non-dyadic results require a separate reference strategy from M2's exact dyadic
arithmetic:

- convert the exact Rational input to a high-guard-precision MPFR computation;
- round once into the target MultiFloat width;
- independently cross-check the pack/result at much higher MPFR precision so the
  candidate is not tested against itself;
- make zero/sign/overflow/subnormal/underflow semantics explicit.

For the arithmetic candidate, use precision-doubling/Newton correction only after
the reference oracle is frozen. Direct FMA/submul residuals should be preferred
where they avoid the intermediate-rounding loss already demonstrated in M5.
Performance is not an M6 acceptance criterion while the safe x5-x8 mul/FMA paths
remain intentionally slow.

## M7 — native linear algebra after competitive arithmetic survivors

Build MultiFloatVec-native DOT/SYRK/GEMM, then TRSM/Cholesky/factor-solve-refine
only after the high-limb arithmetic hot paths are competitive enough that matrix
benchmarks are meaningful. Final acceptance requires downstream solver A/B on
iterations, residuals, certificates, wall time, and memory rather than isolated
microbenchmarks alone.

Every discovered counterexample becomes a permanent regression test.
