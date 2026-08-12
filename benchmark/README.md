# Benchmarks

The benchmark scripts are diagnostic tools, not portable hard gates.

## Scripts

- `smoke.jl`: scalar and Vec4 `fma_fast` versus upstream `x*y+c`;
- `simd_widths.jl`: Vec2/Vec4/Vec8 throughput and speedup;
- `fma3_repair.jl`: exact-value and normalization checks plus direct timing of
  the pre-repair and proof-repaired x3 end networks;
- `division_digits.jl`: rejected quotient-digit candidate versus its full
  residual-product variant and upstream division;
- `codegen.jl`: relative native instruction, stack-reference, call, and FMA
  counts.

Run them from the package root:

```julia
julia --project=. benchmark/smoke.jl
julia --project=. benchmark/simd_widths.jl
julia --project=. benchmark/fma3_repair.jl
julia --project=. benchmark/division_digits.jl
julia --project=. benchmark/codegen.jl
```

## Interpretation

For the FMA scripts, `speedup` or `mul+add/fused` greater than one favors the
candidate. In `fma3_repair.jl`, `repair/old` greater than one is the fixed cost
of the normalization repair. For division, `upstream/specialized` less than one
means upstream is faster.

The scripts use seeded inputs, warm-up calls, repeated samples, and the minimum
observed elapsed time to reduce obvious noise. They do not replace a dedicated
benchmark harness with CPU pinning, frequency control, allocation checks,
confidence intervals, and multiple machines.

A kernel is accepted only after its benchmark result agrees with the intended
workload and an end-to-end downstream A/B. Current decisions are recorded in
`../STATUS.md`.
