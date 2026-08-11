# Verified Float64x5–Float64x8 arithmetic research plan

Research-only until formal verification and numerical gates pass.

## Goal

Develop fixed-length, branch-free arithmetic for 5–8 binary64-component expansions, then use those primitives in a MultiFloat-native linear algebra backend.

Target operations:

- addition/subtraction;
- commutative multiplication;
- fused multiply-add/subtract;
- reciprocal/division/square root;
- SIMD-friendly DOT/SYRK/TRSM/GEMM kernels.

## Acceptance contract

A production candidate must eventually satisfy:

- normalized/non-overlapping output expansion;
- explicit worst-case discarded-tail bound;
- commutativity where mathematically required;
- branch-free hot arithmetic network;
- reproducibility under supported IEEE-754 round-to-nearest semantics;
- adversarial and MPFR/BigFloat differential tests;
- formal FPAN verification when the verifier is capable of proving the required bounds.

## Search strategy

Grow incrementally rather than searching N=8 from scratch:

`N=4 proven seed -> N=5 -> N=6 -> N=7 -> N=8`.

For each N:

1. construct a deliberately safe/over-complete candidate;
2. minimize its gate count and dependency depth;
3. reject quickly using a permanent counterexample corpus;
4. formally verify promising survivors;
5. keep a Pareto frontier in gate count, critical depth, register pressure, proof margin, and measured cycles.

## Counterexample corpus

Include at least:

- leading-limb cancellation;
- multi-limb cancellation;
- values near powers of two and midpoint/tie cases;
- subnormal-adjacent values where supported;
- alternating-sign dot products;
- near-orthogonal Gram/SYRK columns;
- cases extracted from future solver residual/KKT traces.

## Why direct FMA is high priority

Linear algebra hot loops are dominated by `acc += a*b` / `acc -= a*b`. A direct expansion FMA can avoid fully normalizing a multiplication result only to feed it immediately into a second addition/renormalization network.

After `add_N` and `mul_N`, prioritize:

`fma_N/submul_N -> DOT -> SYRK/GEMM -> TRSM/Cholesky`.

Keep a stronger result-relative path for cancellation-sensitive residual/certification work unless the fused network has an adequate result-relative bound.

## Division and square root

Once verified high-limb add/mul exist, prefer precision doubling for Newton/Karp–Markstein-style iterations:

`1 -> 2 -> 4 -> 8`.

Investigate half-to-full and residual-multiply primitives so the final correction step does not pay full x8 cost unnecessarily.

## Milestones

- M0: freeze current x2/x3/x4 behavior and benchmark corpus.
- M1: validate x2/x3/x4 fused FMA + SIMD lane equivalence.
- M2: build a correctness-first x5/x8 reference arithmetic path.
- M3: verified `add5` through `add8`.
- M4: verified commutative `mul5` through `mul8`.
- M5: direct `fma5` through `fma8`.
- M6: `inv/div/sqrt` via precision doubling.
- M7: MultiFloatVec-native DOT/SYRK/TRSM/GEMM.
