import Oscar.Singular as S

const GROBCOV_LIB = "grobcov.lib"

"""
    grobcov(I::Oscar.MPolyIdeal, params, main_vars; options=Int[])

Compute the Groebner cover of the parametric ideal `I` using Singular's `grobcov.lib`.

`params` and `main_vars` are vectors of generators of `base_ring(I)` that specify
which variables are parameters and which are the main polynomial variables.
The polynomial ring is treated as `Q(params)[main_vars]` internally.

Returns a `Vector{GrobCovSegment}`, one element per segment of the partition of
the parameter space.
"""
function grobcov(I::Oscar.MPolyIdeal, params, main_vars; options::Vector{Int} = Int[])
    R_oscar = base_ring(I)
    SF, SR, sing_I = oscar_ideal_to_singular_parametric(I, params, main_vars)

    args = isempty(options) ? Any[sing_I] : Any[sing_I, options]
    raw = S.low_level_caller_ring(GROBCOV_LIB, "grobcov", SR, args)
    # raw is a Vector{Any}; each element is a 3-element vector [lpp, basis, seg]
    _unwrap_and_check(raw, "grobcov")

    return parse_grobcov_output(raw, R_oscar, SF, SR, params, main_vars)
end

"""
    cgsdr(I::Oscar.MPolyIdeal, params, main_vars; options=Int[])

Compute a comprehensive Groebner system with disjoint reduced representation
using Singular's `cgsdr` function from `grobcov.lib`.

Returns a `Vector{CGSSegment}`.
"""
function cgsdr(I::Oscar.MPolyIdeal, params, main_vars; options::Vector{Int} = Int[])
    R_oscar = base_ring(I)
    SF, SR, sing_I = oscar_ideal_to_singular_parametric(I, params, main_vars)

    args = isempty(options) ? Any[sing_I] : Any[sing_I, options]
    raw = S.low_level_caller_ring(GROBCOV_LIB, "cgsdr", SR, args)
    _unwrap_and_check(raw, "cgsdr")

    return parse_cgsdr_output(raw, R_oscar, SF, SR, params, main_vars)
end

# Validate that raw output looks like a list of triples.
function _unwrap_and_check(raw, fname)
    raw isa Vector || error("$fname: unexpected return type $(typeof(raw))")
    isempty(raw) && return
    first_elem = raw[1]
    first_elem isa Vector ||
        error("$fname: expected Vector of triples, got inner type $(typeof(first_elem))")
    length(first_elem) == 3 ||
        error("$fname: expected triples of length 3, got $(length(first_elem))")
end
