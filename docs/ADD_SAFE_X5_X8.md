# Safe Float64x5-x8 addition family

Status: **M3 Experimental correctness baseline**. These functions are not yet
fixed-cost performance kernels.

## Common construction

`Experimental.add_safe(x, y)` supports normalized finite Float64 expansions
with 5 through 8 limbs. Width-specific wrappers `add5_safe` through `add8_safe`
use the same implementation.

For width N the algorithm:

1. applies general `two_sum` to each same-order limb pair;
2. flattens the N exact TwoSum results into an exact 2N-term expansion;
3. fully `renormalize`s the 2N terms;
4. retains the leading N terms;
5. renormalizes that N-limb head.

No FastTwoSum is used, so the safe family has no hidden magnitude-ordering
precondition.

## Exact structural gate

The low N terms of the normalized 2N-term expansion are retained by diagnostics.
For every permanent regression case CI requires the exact dyadic identity

```text
value(result) + value(discarded_tail) = value(x) + value(y).
```

The returned value must also be bitwise identical to the independent M2
`reference_add` oracle, normalized, and bitwise commutative.

The corpus covers ordinary signed random inputs, wide exponent separation,
exact zero/negation identities, deep near-cancellation, powers of two and carry
boundaries, and strongly unbalanced operands.

## Empirical discarded-tail gates

The diagnostic measures

```text
|error| / (u^N * (|x| + |y|)),  u = 2^-53.
```

First accepted hosted-runner measurements:

| Width | Max observed constant | Ordinary nonzero tail | Wide nonzero tail | Near-cancel tail |
|---|---:|---:|---:|---:|
| x5 | 0.0528336 | 583/1000 | 968/1000 | 0/500 |
| x6 | 0.0254805 | 369/600 | 490/600 | 0/200 |
| x7 | 0.0109878 | 335/600 | 401/600 | 0/200 |
| x8 | 0.00541727 | 357/600 | 326/600 | 0/200 |

CI uses **C = 1** independently at every width as a conservative empirical
regression gate. This is not a formal worst-case theorem. The width-specific
measurements are evidence that the safe construction is stable enough to serve
as the reference implementation for the next proof/minimization stage.

## First scalar timing snapshot

Zen 3 hosted runner, Julia 1.10.11, 5,000 additions:

| Width | Safe addition time |
|---|---:|
| x5 | 0.454 ms |
| x6 | 0.545 ms |
| x7 | 0.710 ms |
| x8 | 1.438 ms |

These timings are informational. Performance is not an M3 acceptance criterion.

## Normalized-input contract

The safe family requires normalized inputs. During x6-x8 testing, upstream unary
minus exposed a signed-zero subtlety: direct limb negation can turn trailing
`+0.0` into `-0.0`, while `MultiFloats.isnormalized` uses exact tuple equality.
Identity tests therefore explicitly renormalize negated operands before calling
the strict safe-add API. This does not change the safe-add algorithm or relax its
input contract.

## Next proof target

Before searching for short fixed-cost add5-add8 networks, derive a reviewable
bound on the discarded normalized tail of the common 2N-term construction. Any
future FastTwoSum substitution must prove its source-specific magnitude ordering
and preserve exact-oracle agreement on the permanent corpus.
