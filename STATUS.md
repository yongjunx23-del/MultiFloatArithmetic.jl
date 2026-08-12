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

Two informational CI measurements are now recorded: run 21 used an Ice Lake
runner, while PR run 24 used a Zen 3 runner. The runner change is useful evidence
that these timings are architecture-sensitive rather than portable guarantees.
Full output is preserved in [benchmark/RESULTS.md](benchmark/RESULTS.md).

| Component | Ice Lake evidence | Zen 3 evidence | Decision |
|---|---:|---:|---|
| Float64x2 scalar `fma_fast` | 1.087x | 1.020x | marginal top-level opt-in candidate; no auto-selection |
| Float64x3 scalar `fma_fast` | 0.831x | 0.488x | keep opt-in; scalar regression |
| Float64x4 scalar `fma_fast` | 0.900x | 0.510x | keep opt-in; scalar regression |
| Vec2/Vec4/Vec8 `fma_fast` | 1.033x-2.017x | 0.942x-2.537x | main optimization path, but architecture/width gated |
| specialized expansion × one-limb multiply | 0.981x-1.042x versus full product | 1.007x-1.096x | near-parity experimental helper only |
| quotient-digit division | upstream faster in all six cases | scalar loses badly; Vec4 mixed, including one x3 win | rejected as default; retain for reproducibility |
| Float64x5-Float64x8 arithmetic | none | none | not implemented |
| DOT/SYRK/TRSM/GEMM backend | none | none | not implemented |

For FMA, values above one favor `fma_fast`. For the specialized one-limb
product, values above one favor the specialized product over the zero-padded
full product. Quotient-division evidence is summarized rather than collapsed
into one ratio because its winner changes by width and architecture.

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
