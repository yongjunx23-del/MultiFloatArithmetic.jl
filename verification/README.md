# FPAN verification

This directory mirrors the arithmetic networks in
`src/MultiFloatArithmetic.jl` using the `.fpan` language from
[`dzhang314/FPANVerifier`](https://github.com/dzhang314/FPANVerifier).

## Pinned verifier

CI checks out FPANVerifier commit
`0a1314fa78aeab35793b7354eb097061380982e5` and installs
`z3-solver==4.13.4.0`. Pinning both prevents an upstream lemma, parser, or solver
change from silently changing this repository's proof status.

FPANVerifier expects both Python Z3 bindings and a `z3` executable. CI supplies
the executable through `verification/z3_cli.py`, a minimal adapter that solves
the generated SMT-LIB files with the same pinned Python Z3 package. This avoids
runner-specific APT repositories and eliminates Python/executable version drift.

The upstream verifier currently warns that its SELTZO lemma system is being
redesigned and that precise error-bound/non-overlap proving strength may
temporarily regress. A failed proof at this pin is therefore evidence that the
current tool did not establish the claim; it is not automatically a concrete
floating-point counterexample.

## What is modeled

- normalized x, y, and c input expansions;
- every `two_prod` used by the Julia kernels;
- every rounded one-product term, represented by `two_prod` with its exact
  residual deliberately left unused;
- every rounded `+`, represented by `two_sum` with the error output deliberately
  left unused;
- every explicit `two_sum` and `fast_two_sum` in source order.

Using `two_sum` to model a plain addition is exact for the retained first output:
it is the same rounded sum computed by Julia, while the second output only makes
the discarded rounding residual available to the verifier.

## Current proof obligations

The strict CI command checks:

1. every explicit `fast_two_sum` operand pair satisfies its magnitude
   precondition;
2. the final x2, x3, and x4 output limbs form the expected
   `strongly_dominates` chain.

FPANVerifier validates explicit `fast_two_sum` commands unconditionally. CI does
not pass the verifier's optional `--check-fast-two-sum` flag because that flag
also launches an optimization query for every ordinary `two_sum`, asking whether
it could be replaced by `fast_two_sum`. Those replacement queries do not
strengthen the obligations above and caused the larger source-mirrored networks
to exhaust the solver budget.

Each limb count runs as an independent matrix job with a ten-minute process
limit, so a hard network cannot hide the status of the other two.

These obligations establish structural validity of the accumulation network;
they do **not** yet prove the empirical constants `C_2 = 34`, `C_3 = 184`, and
`C_4 = 812` from `docs/NUMERICAL_CONTRACT.md`, nor do they prove behavior for
NaN, infinity, overflow, or underflow-adjacent inputs.

## Local reproduction

Install `z3-solver==4.13.4.0`, place a compatible `z3` executable on `PATH`,
check out the pinned verifier, then run one or more specs:

```bash
bash verification/run_fpan.sh /path/to/FPANVerifier verification/fma2.fpan
```

Any line beginning with `ERROR:` makes the wrapper fail, because FPANVerifier
currently reports refuted `prove` statements and invalid `fast_two_sum` uses in
text while continuing to process the remaining file.
