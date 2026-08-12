# Verified Float64x5–Float64x8 arithmetic research plan

Research-only until formal verification and numerical gates pass.

## Goal

Develop reliable fixed-length arithmetic for 5–8 binary64-component expansions,
then optimize accepted primitives for a MultiFloat-native linear algebra backend.

Target operations:

- addition/subtraction;
- commutative multiplication;
- fused multiply-add/subtract;
- reciprocal/division/square root;
- SIMD-friendly DOT/SYRK/TRSM/GEMM kernels.

## Acceptance contract

A production performance candidate must eventually satisfy:

- normalized/non-overlapping output expansion;
- explicit worst-case discarded-tail bound;
- commutativity where mathematically required;
- reproducibility under supported IEEE-754 round-to-nearest semantics;
- adversarial and exact-oracle differential tests;
- all FastTwoSum magnitude assumptions proved rather than inferred from typical
  inputs;
- formal verification when the available verifier can establish the required
  properties;
- branch-free/fixed-cost execution only after the correctness network is sound.

The x4 cancellation defect permanently changed the search order: correctness and
normalization are established before branch-free minimization.

## M2 exact oracle — complete baseline

`MultiFloatArithmetic.Experimental` provides `reference_add`, `reference_sub`,
`reference_mul`, and `reference_fma` for Float32/Float64 x5-x8.

Each scalar input is converted to exact `Rational{BigInt}` dyadics, the operation
is performed exactly, and the result is packed once. CI cross-checks the internal
packing precision against independent 8192-bit BigFloat packing. Vector reference
methods are lane-wise scalar calls.

This oracle is the rejection standard for M3-M5 and is not a performance target.

## Search strategy

Grow incrementally:

`safe x5 -> validated x5 -> safe x6 -> safe x7 -> safe x8 -> fixed-cost search`.

For each N and operation:

1. construct a deliberately safe/over-complete candidate;
2. compare against the exact M2 oracle on deterministic/adversarial inputs;
3. measure discarded-tail/error constants, normalization failures, symmetry, and
   lane equivalence;
4. reject immediately on any structural counterexample;
5. establish source-specific ordering bounds before introducing FastTwoSum;
6. formally verify promising survivors where feasible;
7. only then minimize gate count/dependency depth and benchmark;
8. keep a Pareto frontier in gate count, critical depth, register pressure,
   proof margin, and measured cycles.

## Counterexample corpus

Include at least:

- leading- and multi-limb cancellation;
- exact product/add cancellation;
- powers of two, carry boundaries, and midpoint/tie neighborhoods;
- strongly unbalanced operand magnitudes;
- subnormal-adjacent values where supported;
- alternating-sign dot products;
- near-orthogonal Gram/SYRK columns;
- cases extracted from future solver residual/KKT traces.

Every discovered counterexample becomes a permanent regression test.

## M3 addition strategy

### add5_safe — accepted Experimental baseline

The first Float64x5 candidate uses no FastTwoSum:

1. five same-index `two_sum` operations;
2. full renormalization of the exact ten-term expansion;
3. retain the leading five terms and renormalize the head.

The current permanent corpus requires:

- exact identity between returned head + discarded normalized tail and `x+y`;
- bitwise equality with `reference_add`;
- normalized output;
- bitwise commutativity;
- exact identities, deep cancellation, power-of-two/carry boundaries, and
  strongly unbalanced operands.

The first 2,500-case diagnostic found zero oracle/normalization/commutativity
failures and a maximum observed
`|err|/(u^5(|x|+|y|)) ≈ 0.05284`. CI uses `C=1` as a conservative empirical
regression gate, not a proof.

### Next

Generalize the same construction to safe add6/add7/add8 and measure a separate
empirical tail constant for each width. Do not minimize add5 yet: first establish
that the construction scales cleanly across all M3 widths.

After safe add5-add8 exist, derive a formal discarded-tail bound and only then
search fixed-cost compression networks.

## Why direct FMA is high priority

Linear algebra hot loops are dominated by `acc += a*b` / `acc -= a*b`. A direct
expansion FMA can avoid fully normalizing a multiplication result only to feed it
immediately into a second addition/renormalization network.

After `add_N` and `mul_N`, prioritize:

`fma_N/submul_N -> DOT -> SYRK/GEMM -> TRSM/Cholesky`.

Keep a stronger result-relative path for cancellation-sensitive residual/
certification work unless the fused network has an adequate result-relative
bound.

## Division and square root

Once verified high-limb add/mul exist, prefer precision doubling for
Newton/Karp–Markstein-style iterations:

`1 -> 2 -> 4 -> 8`.

Investigate half-to-full and residual-multiply primitives so the final correction
step does not pay full x8 cost unnecessarily.

## Milestones

- M0: **complete** — initial x2/x3/x4 behavior and benchmark corpus.
- M1/M1.1: **complete baseline** — x2/x3 structural verification, x3 repair,
  cancellation-safe x4 baseline.
- M2: **complete baseline** — exact x5-x8 reference add/sub/mul/FMA oracle.
- M3: **in progress** — `add5_safe` accepted Experimental baseline; add6-add8
  safe baselines next, followed by formal tail bound and fixed-cost search.
- M4: verified commutative `mul5` through `mul8`.
- M5: direct `fma5` through `fma8`.
- M6: `inv/div/sqrt` via precision doubling.
- M7: MultiFloatVec-native DOT/SYRK/TRSM/GEMM.
