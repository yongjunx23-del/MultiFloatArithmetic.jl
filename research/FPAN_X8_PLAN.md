# Verified Float64x5–Float64x8 arithmetic research plan

Research-only until formal verification and numerical gates pass.

## Goal

Develop reliable fixed-length arithmetic for 5–8 binary64-component expansions,
then optimize the accepted primitives for a MultiFloat-native linear algebra
backend.

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
- formal FPAN verification when the verifier is capable of proving the required
  bounds;
- branch-free/fixed-cost execution only after the correctness network is sound.

The x4 cancellation defect changed the search order deliberately: correctness
and normalization are now established before branch-free minimization.

## M2 exact oracle — implemented

The package now has an independent x5-x8 oracle under
`MultiFloatArithmetic.Experimental`:

- `reference_add`;
- `reference_sub`;
- `reference_mul`;
- `reference_fma`.

It converts each finite normalized MultiFloat input to exact `Rational{BigInt}`
dyadics, performs the operation exactly, and packs the result once. The internal
packing precision is cross-checked against 8192-bit BigFloat packing in CI.
`MultiFloatVec` reference methods are lane-wise scalar calls, intentionally
avoiding a shared vector arithmetic implementation.

This oracle is the rejection standard for M3-M5 and is not a performance target.

## Search strategy

Grow incrementally rather than searching N=8 from scratch:

`safe x5 -> validated x5 -> x6 -> x7 -> x8`.

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

- leading-limb cancellation;
- multi-limb cancellation;
- exact product/add cancellation;
- values near powers of two and midpoint/tie cases;
- subnormal-adjacent values where supported;
- alternating-sign dot products;
- near-orthogonal Gram/SYRK columns;
- cases extracted from future solver residual/KKT traces.

Every discovered counterexample becomes a permanent regression test.

## M3 addition strategy

Start with x5 only. The first candidate should prioritize an explicit
renormalized expansion and measurable discarded tail rather than a hand-minimized
FastTwoSum network.

Acceptance sequence:

`add5_safe -> adversarial differential test -> tail-bound study -> structural proof -> add5_fast -> add6...add8`.

Do not promote a short addition network merely because ordinary random inputs
match the oracle.

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

- M0: **complete** — freeze initial x2/x3/x4 behavior and benchmark corpus.
- M1/M1.1: **complete baseline** — validate x2/x3; repair x3; replace unsafe x4
  final compression with the cancellation-safe baseline.
- M2: **implemented, validation/merge in progress** — exact x5-x8 reference
  add/sub/mul/FMA oracle.
- M3: next — verified `add5` through `add8`, beginning with safe x5 addition.
- M4: verified commutative `mul5` through `mul8`.
- M5: direct `fma5` through `fma8`.
- M6: `inv/div/sqrt` via precision doubling.
- M7: MultiFloatVec-native DOT/SYRK/TRSM/GEMM.
