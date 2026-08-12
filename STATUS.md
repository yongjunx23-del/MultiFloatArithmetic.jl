# Project status

Status date: **2026-08-12**
Current milestone: **M1 — x2/x3/x4 FMA empirical validation and API freeze**

## Executive assessment

The repository is a coherent research prototype with passing correctness CI,
SIMD lane-equivalence tests, adversarial cancellation coverage, code-generation
diagnostics, and benchmark scripts. It is not yet a production arithmetic
library because the stated formal verification gate has not been completed and
most of the x5-x8/linear-algebra roadmap is still unimplemented.

A practical estimate is:

- M0/M1 engineering and empirical validation: substantially complete;
- formal proof package for current kernels: pending;
- M2-M6 higher-limb arithmetic: not started;
- M7 MultiFloat-native linear algebra: not started.

## Kernel decisions

| Component | Evidence | 2026-08-11 performance evidence | Decision |
|---|---|---:|---|
| Float64x2 scalar `fma_fast` | BigFloat differential tests, normalization, commutativity, cancellation corpus | 1.087x vs `x*y+c` | top-level opt-in candidate |
| Float64x3 scalar `fma_fast` | same | 0.831x | keep opt-in; do not auto-select |
| Float64x4 scalar `fma_fast` | same | 0.900x | keep opt-in; do not auto-select |
| Vec2/Vec4/Vec8 `fma_fast` | scalar-lane bitwise equivalence and operand-relative checks | 1.033x-2.017x across tested cases | continue as the main optimization path |
| specialized expansion × one-limb multiply | bitwise equivalence to zero-padded upstream `mfmul` | only about 0.6%-4.2% over the full residual product | experimental helper only |
| quotient-digit division | BigFloat relative-error checks and lane equivalence | upstream `/` was about 1.09x-4.22x faster | rejected as default; retain for reproducibility |
| Float64x5-Float64x8 arithmetic | research plan only | none | not implemented |
| DOT/SYRK/TRSM/GEMM backend | roadmap only | none | not implemented |

The timing data above came from CI run 21 on an Ice Lake hosted runner. These
numbers are directional, not portable performance guarantees.

## Acceptance gates

A kernel may move from `Experimental` to the top-level API only when all
applicable gates pass:

1. **Domain and semantics:** finite-input domain, rounding mode, underflow
   assumptions, normalization invariant, and error metric are explicit.
2. **Deterministic correctness:** permanent adversarial cases, seeded randomized
   tests, BigFloat/MPFR differential checks, and SIMD lane equivalence pass.
3. **Verification:** non-overlap and discarded-tail/error bounds have a
   reviewable proof or machine-verifier artifact.
4. **Code generation:** no hidden calls, accidental reassociation, or unexpected
   stack/register explosion on supported targets.
5. **Performance:** repeated architecture-specific measurements show a stable
   benefit for the intended scalar or vector workload.
6. **Downstream value:** an end-to-end solver A/B improves time or memory without
   degrading iterations, residuals, or certificates.

## Immediate next milestone

The next technically meaningful milestone is **M1.1: formalize and verify the
x2/x3/x4 accumulation networks**, while preserving the current Julia code as a
frozen executable specification. Only after that should the project start a
correctness-first x5 reference implementation.
