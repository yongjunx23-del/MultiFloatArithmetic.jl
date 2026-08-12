# Compatibility shim for the short-lived M6 `reference_recip` spelling.
#
# `reference_inverse_div_sqrt.jl` is the sole authoritative non-dyadic oracle:
# it adaptively doubles MPFR precision until the final MultiFloat pack stabilizes
# at two consecutive precisions and also provides division and square root.
#
# Keep only the reciprocal spelling alias here so the already-created x5 Newton
# candidate and any Experimental callers do not break. Do NOT define
# `reference_div` here; having two implementations caused method-overwrite errors
# during Julia module precompilation.

"""
    reference_recip(x)

Compatibility alias for [`reference_inv`](@ref). New M6 work should use
`reference_inv`; both spellings execute the same adaptive stabilized oracle.
"""
reference_recip(x::MultiFloat{T,N}) where {T,N} = reference_inv(x)

function reference_recip(x::MultiFloatVec{W,T,N}) where {W,T,N}
    return reference_inv(x)
end
