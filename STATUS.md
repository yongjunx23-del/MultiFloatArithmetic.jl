# Project status

Status date: **2026-08-12**
Current milestone: **M1.1 — x2/x3 structural verification and x4 decomposition**

## Executive assessment

The repository is now a reviewable v0.1 research package rather than a loose
prototype: the public/experimental API boundary is explicit, conventional Julia
CI is cross-platform and bounds-checked, negative performance results are
recorded, and source-mirrored formal-verification assets exist. It is still not
a production arithmetic library because the empirical error constants are not
machine-proved, x4 structural verification is unresolved, and the x5-x8 and
linear-algebra roadmap remains unimplemented.

A practical milestone assessment is:

- M0/M1 engineering, API freeze, and empirical validation: substantially
  complete;
- x2 structural proof: complete at the pinned verifier/toolchain;
- x3 structural proof: complete after the proof-driven final normalization
  repair;
- x4 structural proof: unresolved because the universal model exceeds the
  hosted-runner budget, with no refutation observed;
- formal proof of `C_2`, `C_3`, and `C_4`: pending;
- M2-M6 higher-limb arithmetic: not started;
- M7 MultiFloat-native linear algebra: not started.

## Formal finding and repair

The original x3 end network proved `z1 strongly_dominates z2` but FPANVerifier
refuted `z0 strongly_dominates z1` in an abstract precision-11 model. Randomized
Julia tests did not expose a concrete non-normalized output, so the model was
treated as a structural gap rather than presented as a concrete IEEE input.

A minimal A/B compared three fixed-cost repairs. A single final `two_sum(z0,
z1)` was insufficient because the resulting tail relation was refuted. Adding
`fast_two_sum(z1, z2)` after that operation proved both final relations. A full
fixed renormalization pass also proved them but uses a more general operation on
the tail. The selected repair is therefore the two-operation candidate.

## Kernel decisions

Two informational CI measurements are recorded: run 21 used an Ice Lake runner,
while PR run 24 used a Zen 3 runner. The runner change demonstrates that these
timings are architecture-sensitive rather than portable guarantees. Full output
is preserved in [benchmark/RESULTS.md](benchmark/RESULTS.md).

| Component | Ice Lake evidence | Zen 3 evidence | Decision |
|---|---:|---:|---|
| Float64x2 scalar `fma_fast` | 1.087x | 1.020x | marginal top-level opt-in candidate; no auto-selection |
| Float64x3 scalar `fma_fast` | 0.831x | 0.488x | keep opt-in; scalar regression; repair overhead measured separately |
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

A kernel may be described as fully verified only when all applicable gates pass:

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

Finish **M1.1** by decomposing or fixed-precision validating the x4 network and
by deriving machine-checkable error constants for x2/x3. Only after those tasks
should the project start a correctness-first x5 reference implementation.
