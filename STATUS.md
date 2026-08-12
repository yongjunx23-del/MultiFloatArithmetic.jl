# Project status

Status date: **2026-08-12**
Current milestone: **M1.1 — x2/x3 structural verification and x4 correctness hardening**

## Executive assessment

The repository is now a reviewable v0.1 research package with an explicit
public/experimental API boundary, cross-platform bounds-checked CI, adversarial
BigFloat differential tests, SIMD lane-equivalence tests, benchmark records, and
pinned FPANVerifier assets.

It is **not** a production arithmetic library yet. The current state is:

- x2 FMA structural proof: complete at the pinned verifier/toolchain;
- x3 FMA structural proof: complete after the proof-driven final-normalization
  repair;
- x4 FMA: original fixed-cost end network has a confirmed cancellation bug;
  a conservative `renormalize` baseline now removes the invalid FastTwoSum
  assumptions and passes the empirical correctness suite;
- formal proof of the empirical constants `C_2`, `C_3`, and `C_4`: pending;
- x5-x8 arithmetic: not implemented;
- MultiFloat-native DOT/SYRK/TRSM/GEMM backend: not implemented.

## x4 finding

The original x4 final compression began with
`fast_two_sum(b, a1)`. FPAN decomposition refuted the required magnitude
ordering under cancellation. A concrete Float64 diagnostic then made the issue
observable in ordinary Julia execution.

For 10,000 seeded cases with `c` chosen as the x4 representation of `-x*y`:

- first-FastTwoSum ordering violations: **3,755 / 10,000**;
- old x4 non-normalized outputs: **5,258 / 10,000**;
- old output differed from the safe baseline: **5,258 / 10,000**;
- safe-baseline non-normalized outputs: **0 / 10,000**.

On 20,000 ordinary random values the old and safe x4 outputs were bitwise
identical, which explains why the earlier broad randomized suite did not expose
the cancellation defect.

## Current x4 baseline

The correctness-first x4 path now:

1. uses general `two_sum` instead of the two risky final-compression
   `fast_two_sum` operations;
2. finishes with `MultiFloats.renormalize`;
3. preserves the existing x4 product/accumulation network before final
   compression.

This intentionally gives up the package's fixed-cost branch-free property for
x4 until a source-specific fixed-cost compression is proved. Brute-force
assumption-free fixed-tail A/B candidates were rejected: at p=24 and p=53 they
still failed the top two non-overlap obligations.

## Performance snapshot

GitHub-hosted runner, Julia 1.10.11, `Sys.CPU_NAME = generic`:

| Workload | Old x4 | Safe x4 | Relative result |
|---|---:|---:|---:|
| scalar, 20k | 0.062 ms | 0.345 ms | safe/old 5.593x |
| Vec2, 20k lanes | 0.155 ms | 0.211 ms | safe/old 1.364x |
| Vec4, 20k lanes | 0.080 ms | 0.109 ms | safe/old 1.352x |
| Vec8, 20k lanes | 0.064 ms | 0.068 ms | safe/old 1.056x |

Against upstream `x*y+c`, the safe x4 SIMD smoke still measured about **1.59x
(Vec2), 1.59x (Vec4), and 1.31x (Vec8)** on this runner. Scalar x4 is not a
performance candidate.

## Kernel decisions

| Component | Decision |
|---|---|
| x2 scalar `fma_fast` | marginal architecture-dependent opt-in candidate |
| x3 scalar `fma_fast` | correctness-verified but slower; do not auto-select |
| x2/x3 SIMD `fma_fast` | continue architecture-specific optimization path |
| x4 scalar `fma_fast` | correctness baseline only; upstream is faster |
| x4 SIMD `fma_fast` | safe baseline retains useful speedup; fixed-cost repair still research |
| quotient-digit division | rejected as default; retain under `Experimental` |
| expansion × one-limb multiply | experimental helper only |
| x5-x8 | not implemented |
| native linear algebra backend | not implemented |

## Acceptance gates

A kernel may be described as fully verified only when all applicable gates pass:

1. domain, rounding, underflow, normalization, and error metric are explicit;
2. deterministic/adversarial tests and BigFloat differential checks pass;
3. structural non-overlap and FastTwoSum preconditions are proved or avoided;
4. discarded-tail/error constants have reviewable proof evidence;
5. code generation is inspected on supported targets;
6. performance benefit is repeatable on the intended architecture/workload;
7. downstream solver A/B improves time or memory without degrading residuals or
   certificates.

## Immediate next milestone

Finish M1.1 by making the x4 safe baseline part of the frozen source contract,
then derive machine-checkable error bounds for x2/x3 and a source-specific x4
compression proof. Do **not** start x5 optimization by copying the old x4
FastTwoSum assumptions.
