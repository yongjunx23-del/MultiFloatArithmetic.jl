# Float64x5 reciprocal correctness baseline

M6 starts with a correctness-first reciprocal for `MultiFloat{Float64,5}`.
The accepted Experimental entrypoint is `recip5_safe`.

## Oracle

The authoritative non-dyadic oracle is `reference_inv`. It keeps the input exact
as `Rational{BigInt}`, evaluates the reciprocal with MPFR, and doubles MPFR
precision until the final five-limb pack is identical at two consecutive
precisions (up to 32768 bits). Permanent tests independently validate the oracle
against a 16384-bit pack on a separate corpus.

`reference_recip` remains only a compatibility spelling for `reference_inv`.
There is no second reciprocal implementation behind that alias.

## Accepted construction

The accepted x5 path is precision-doubling rather than repeated x5 correction:

1. truncate the normalized x5 input to its leading x4 expansion;
2. compute the reciprocal with upstream MultiFloats.jl x4 `inv`;
3. lift the normalized x4 reciprocal to x5 by adding a zero fifth limb;
4. compute the x5 residual directly as `x*z - 1` with `fma5_safe`;
5. perform one Newton correction directly as `z - z*(x*z - 1)` with a second
   `fma5_safe`.

The direct-FMA correction is intentional. M5 showed that rounded multiplication
followed by rounded addition is not a valid substitute for direct FMA under
cancellation.

## A/B against a Float64 seed

The rejected slower alternative starts from `inv(x._limbs[1])` and applies three
x5 direct-FMA Newton corrections. Both paths are retained in the benchmark as an
A/B, but only the x4-seeded one-correction path backs `recip5_safe`.

First accepted Zen 3 / Julia 1.10.11 snapshot, 120 seeded inputs:

| Path | Oracle mismatches | Worst exact residual bits | Time / 120 |
|---|---:|---:|---:|
| x4 seed + 1 x5 correction | 0 | 269.1 | 19.450 ms |
| Float64 seed + 3 x5 corrections | 0 | 269.1 | 57.494 ms |

The x4 seed itself already had about 209.8 worst-case exact residual bits on this
corpus. One x5 Newton correction raised that to about 269.1 bits and matched
`reference_inv` bit-for-bit. The three-correction path was about **2.96x slower**
with no precision advantage in the accepted corpus.

These timings use the current over-complete direct-FMA correctness baseline, so
they are not a production reciprocal performance claim.

## Permanent gates

The x5 reciprocal corpus requires:

- bitwise equality with adaptive `reference_inv`;
- normalized finite output;
- exact `Rational{BigInt}` residual contraction;
- end-to-end relative-error telemetry tighter than `2^-250`;
- identity/sign/power-of-two cases;
- explicit zero/nonfinite/unnormalized and seed-domain handling;
- cross-platform Linux Julia 1.10/current and macOS current coverage.

The accepted branch passes all gates.

## Current domain

`recip5_safe` currently requires a finite normalized nonzero x5 input whose
leading Float64 limb is normal and whose leading-limb reciprocal is also normal.
It additionally inherits the conservative product-side underflow exclusions of
`fma5_safe`.

Subnormal-leading inputs and inputs such as `floatmax(Float64)` whose hardware
leading-limb reciprocal is subnormal are deliberately rejected until gradual
underflow is handled explicitly.

## Next

- establish the corresponding x8 direct-FMA reciprocal correction count;
- then build division with a quotient residual `x - q*y` evaluated by direct FMA
  and require bitwise `reference_div` agreement;
- keep reciprocal/division correctness work separate from the fixed-cost x8
  multiplication/FMA performance-reduction track.
