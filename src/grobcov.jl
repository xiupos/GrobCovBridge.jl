import Oscar.Singular as S

const GROBCOV_LIB = "grobcov.lib"

"""
    grobcov(I::Oscar.MPolyIdeal, params, main_vars; can=1, ext=0, comment=0)

Compute the Groebner cover of the parametric ideal `I` using Singular's `grobcov.lib`.

`params` and `main_vars` are vectors of generators of `base_ring(I)` that specify
which variables are parameters and which are the main polynomial variables.
The polynomial ring is treated as `Q(params)[main_vars]` internally.

Options (mirror grobcov.lib named options):
- `can`:     0 or 1 (default 1). `can=1` homogenizes the full ideal before computing
             the Groebner Cover (canonical result); `can=0` only homogenizes the basis.
- `ext`:     0 or 1 (default 0). `ext=1` computes the full representation of each basis
             element (a sheaf of polynomials); `ext=0` gives the generic representation.
- `comment`: 0–3 (default 0). Higher values print progress information from Singular.

Returns a `Vector{GrobCovSegment}`, one element per segment of the partition of
the parameter space.
"""
function grobcov(I::Oscar.MPolyIdeal, params, main_vars;
                 can::Int = 1, ext::Int = 0, comment::Int = 0)
    R_oscar = base_ring(I)
    SF, SR, sing_I = oscar_ideal_to_singular_parametric(I, params, main_vars)

    args = Any[sing_I, "can", can, "ext", ext, "comment", comment]
    raw = S.low_level_caller_ring(GROBCOV_LIB, "grobcov", SR, args)
    _unwrap_and_check(raw, "grobcov")

    return parse_grobcov_output(raw, R_oscar, SF, SR, params, main_vars)
end

"""
    cgsdr(I::Oscar.MPolyIdeal, params, main_vars; can=2, out=1, comment=0)

Compute a comprehensive Groebner system with disjoint reduced representation
using Singular's `cgsdr` function from `grobcov.lib`.

Options:
- `can`:     0, 1, or 2 (default 2). Controls homogenization strategy.
             `can=2` uses the given basis directly (KSW algorithm only, fastest).
             `can=0` homogenizes the basis; `can=1` homogenizes the full ideal.
- `out`:     0 or 1 (default 1). `out=1` returns segments as `(E, N, basis)` triples
             (V(E) \\ V(N)). `out=0` groups by lpp; not compatible with `can=0,1`.
- `comment`: 0 or 1 (default 0). Set to 1 to print progress information.

Returns a `Vector{CGSSegment}`.
"""
function cgsdr(I::Oscar.MPolyIdeal, params, main_vars;
               can::Int = 2, out::Int = 1, comment::Int = 0)
    # grobcov.lib forces out=0 when can<2, which produces a different output structure
    # that the current parser does not handle.
    if can < 2 && out != 0
        throw(ArgumentError(
            "cgsdr: can=$can forces out=0 inside grobcov.lib. " *
            "Pass out=0 explicitly, or use the default can=2 for (E,N,basis) output."))
    end

    R_oscar = base_ring(I)
    SF, SR, sing_I = oscar_ideal_to_singular_parametric(I, params, main_vars)

    args = Any[sing_I, "can", can, "out", out, "comment", comment]
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
