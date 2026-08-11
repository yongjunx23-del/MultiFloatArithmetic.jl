using InteractiveUtils
using MultiFloatArithmetic
using MultiFloats

@inline fused(x, y, c) = fma_fast(x, y, c)
@inline separate(x, y, c) = x * y + c

function native_metrics(f, ::Type{T}) where {T}
    io = IOBuffer()
    code_native(io, f, Tuple{T,T,T}; syntax=:intel, debuginfo=:none)
    asm = String(take!(io))
    lines = split(asm, '\n')

    instructions = String[]
    for line in lines
        s = strip(line)
        isempty(s) && continue
        startswith(s, ".") && continue
        endswith(s, ":") && continue
        startswith(s, ";") && continue
        # Intel-syntax instruction lines emitted by Julia begin with a mnemonic
        # after whitespace. Keep this deliberately approximate: the diagnostic is
        # for comparing the two implementations under the same compiler/runner.
        occursin(r"^[A-Za-z][A-Za-z0-9_.]*[ \t]", s) || continue
        push!(instructions, s)
    end

    stack_refs = count(s -> occursin("rsp", lowercase(s)) || occursin("rbp", lowercase(s)), instructions)
    calls = count(s -> startswith(lowercase(s), "call"), instructions)
    fmas = count(s -> occursin("vfm", lowercase(s)) || occursin("fma", lowercase(s)), instructions)

    return (
        instructions=length(instructions),
        stack_refs=stack_refs,
        calls=calls,
        fma_instructions=fmas,
        asm_bytes=ncodeunits(asm),
    )
end

function show_pair(label, ::Type{T}) where {T}
    mf = native_metrics(fused, T)
    ms = native_metrics(separate, T)
    println(label)
    println("  fused:    ", mf)
    println("  mul+add:  ", ms)
end

println("Codegen diagnostics; relative indicators only")
for (N, T) in enumerate((Float64x2, Float64x3, Float64x4))
    show_pair("$(T) scalar", T)
    V = MultiFloatVec{4,Float64,N}
    show_pair("$(T) Vec4", V)
end
