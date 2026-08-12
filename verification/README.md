# FPAN verification

This directory contains the machine-checkable structural models for the
**fixed-cost x2 and x3** `fma_fast` networks in `src/MultiFloatArithmetic.jl`.
The accepted x4 path is intentionally different: it ends in upstream
`MultiFloats.renormalize`, so x4 normalization is not represented as a fixed
FPAN network.

## Pinned verifier

CI checks out FPANVerifier commit
`0a1314fa78aeab35793b7354eb097061380982e5` and installs
`z3-solver==4.13.4.0`. Pinning both prevents an upstream lemma, parser, or solver
change from silently changing this repository's proof status.

FPANVerifier expects both Python Z3 bindings and a `z3` executable. CI supplies
the executable through `verification/z3_cli.py`, a minimal adapter that solves
the generated SMT-LIB files with the same pinned Python Z3 package.

The upstream verifier warns that its SELTZO lemma system is under redesign. A
failed proof at this pin is evidence that the current model did not establish a
claim; it is not automatically a concrete IEEE counterexample.

## x2

`fma2.fpan` mirrors the x2 source network and proves the explicit FastTwoSum
magnitude precondition plus the final output non-overlap relation.

## x3

The original x3 end network had a structural normalization gap. The selected
repair adds

```text
two_sum(z0, z1)
fast_two_sum(z1, z2)
```

and `fma3.fpan` proves both final non-overlap relations as well as the explicit
FastTwoSum preconditions. Diagnostic A/B history is retained in closed PR #9.

## x4

The historical fixed-cost x4 model exposed an invalid
`fast_two_sum(b, a1)` ordering assumption under destructive cancellation.
Concrete Float64 tests then reproduced non-normalized old outputs. The accepted
x4 source therefore uses general TwoSum transforms and delegates final
normalization to `MultiFloats.renormalize`.

`fma4.fpan` is retained only as a **historical rejected model** so the failure is
auditable; it is not run by CI and must not be described as source-mirrored
verification of the accepted x4 implementation. The x4 decomposition and repair
experiments are retained in closed PR #10 and its Actions artifacts.

## CI policy

Only x2 and x3 are strict FPAN jobs. x4 correctness is guarded by the normal
Julia suite, including BigFloat differential checks, SIMD lane equivalence,
destructive cancellation, and the dedicated x4 safe-baseline diagnostic.

FPANVerifier validates explicit `fast_two_sum` commands unconditionally. CI does
not pass `--check-fast-two-sum` because that optional flag additionally asks
whether every ordinary TwoSum can be optimized to FastTwoSum and makes the model
needlessly expensive.

These structural obligations do **not** prove the empirical constants
`C_2 = 34`, `C_3 = 184`, and `C_4 = 812`, nor NaN/infinity/overflow/
underflow-adjacent semantics.

## Local reproduction

```bash
bash verification/run_fpan.sh /path/to/FPANVerifier verification/fma2.fpan
bash verification/run_fpan.sh /path/to/FPANVerifier verification/fma3.fpan
```

Any line beginning with `ERROR:` makes the wrapper fail.
