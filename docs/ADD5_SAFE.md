# Safe Float64x5 addition baseline

Status: **M3 Experimental correctness baseline**. This is not yet a fixed-cost
performance kernel.

## Algorithm

`Experimental.add5_safe(x, y)` accepts finite normalized `Float64x5` values and:

1. runs `two_sum` on each same-order limb pair, producing ten terms whose exact
   sum equals `x + y`;
2. runs `MultiFloats.renormalize` on all ten terms;
3. retains the leading five normalized terms;
4. renormalizes the five-term head once more.

No `fast_two_sum` is used, so the baseline has no hidden magnitude-ordering
precondition.

## Exact tail accounting

The test path retains the low five terms of the normalized ten-term expansion.
For every permanent regression case CI requires the exact dyadic identity

```text
value(add5_safe(x,y)) + value(discarded_tail) = value(x) + value(y).
```

This makes the truncation error auditable rather than inferred from a floating
reference calculation.

## Oracle gate

The returned five-limb value must also be bitwise identical to the independent
M2 `reference_add` oracle. The current corpus covers:

- ordinary signed random inputs;
- wide exponent separation;
- exact zero/negation identities;
- deep near-cancellation;
- powers of two and carry boundaries;
- strongly unbalanced operand magnitudes.

The corpus currently shows no oracle mismatch, normalization failure, or
commutativity failure.

## Empirical error gate

The diagnostic records

```text
|error| / (u^5 * (|x| + |y|)),  u = 2^-53.
```

The first accepted corpus measured a maximum of approximately **0.05284**.
CI uses **C = 1** only as a conservative empirical regression gate; it is not a
formal theorem. A future proof should replace this empirical constant with a
reviewable bound derived from the normalized ten-term tail.

Near-cancellation cases in the first corpus were represented exactly by the
five-limb result, but no strong result-relative guarantee is claimed from that
observation alone.

## Performance role

The implementation intentionally uses iterative `renormalize`. It is the
correctness baseline for later fixed-cost `add5` work, not the final throughput
path.

Any optimized replacement must preserve:

- bitwise agreement with `reference_add` on the permanent corpus;
- exact identities and commutativity;
- normalized output;
- an explicit discarded-tail/error bound;
- proof of every `fast_two_sum` magnitude assumption it introduces.
