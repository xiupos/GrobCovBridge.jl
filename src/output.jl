"""
    GrobCovSegment

One segment from the output of `grobcov`. The parameter space is partitioned into
segments; on each segment the Groebner cover basis `basis` has leading power product
list `lpp`.

Fields:
- `lpp`:     leading power products of the basis (Oscar MPolyRingElem vector)
- `basis`:   reduced Groebner basis as an Oscar MPolyIdeal
- `segment`: P-representation of the parameter segment as
             `[[prime_ideal, [excl_ideal_1, ...]], ...]`
             meaning ∪_i (V(prime_i) \\ ∪_j V(excl_ij))
"""
struct GrobCovSegment
    lpp::Vector{<:Oscar.MPolyRingElem}
    basis::Oscar.MPolyIdeal
    segment::Vector  # Vector of [MPolyIdeal, Vector{MPolyIdeal}]
end

"""
    CGSSegment

One segment from the output of `cgsdr`. Each segment gives conditions on the
parameters (null and nonnull ideals) and a basis for the ideal on that segment.

Fields:
- `null`:    ideal whose variety is the null set of the segment (Oscar MPolyIdeal)
- `nonnull`: ideal whose variety must be avoided (Oscar MPolyIdeal)
- `basis`:   Groebner basis on this segment (Oscar MPolyIdeal)
"""
struct CGSSegment
    null::Oscar.MPolyIdeal
    nonnull::Oscar.MPolyIdeal
    basis::Oscar.MPolyIdeal
end

# Parse grobcov raw output (Vector{Any} from low_level_caller_ring) into GrobCovSegments.
# raw_result: result[i] = [lpp_ideal, basis_ideal, seg_desc]
#   lpp_ideal:   sideal of leading power products (monomials)
#   basis_ideal: sideal of basis polynomials
#   seg_desc:    Vector{Any} P-representation [[prime, [excl,...]], ...]
# R_oscar:   the Oscar MPolyRing the user passed in
# SF:        Singular N_FField (parameter field)
# SR:        Singular PolyRing (main variables)
# params:    parameter variables in R_oscar
# main_vars: main variables in R_oscar
function parse_grobcov_output(raw_result, R_oscar, SF, SR, params, main_vars)
    result = GrobCovSegment[]
    for triple in raw_result
        lpp_ideal   = triple[1]  # sideal
        basis_ideal = triple[2]  # sideal
        seg_raw     = triple[3]  # Vector{Any}

        lpp_polys = _sideal_to_oscar_vec(lpp_ideal, SF, SR, R_oscar, params, main_vars)
        basis_polys = _sideal_to_oscar_vec(basis_ideal, SF, SR, R_oscar, params, main_vars)
        basis_oscar = ideal(R_oscar, basis_polys)

        seg = _parse_segment(seg_raw, SF, SR, R_oscar, params, main_vars)
        push!(result, GrobCovSegment(lpp_polys, basis_oscar, seg))
    end
    return result
end

# Parse cgsdr raw output into CGSSegments.
# raw_result: result[i] = [E_ideal, N_ideal, B_ideal]
function parse_cgsdr_output(raw_result, R_oscar, SF, SR, params, main_vars)
    result = CGSSegment[]
    for triple in raw_result
        e_ideal = triple[1]  # sideal (null conditions)
        n_ideal = triple[2]  # sideal (nonnull conditions)
        b_ideal = triple[3]  # sideal (basis)

        e_polys = _sideal_to_oscar_vec(e_ideal, SF, SR, R_oscar, params, main_vars)
        n_polys = _sideal_to_oscar_vec(n_ideal, SF, SR, R_oscar, params, main_vars)
        b_polys = _sideal_to_oscar_vec(b_ideal, SF, SR, R_oscar, params, main_vars)

        push!(result, CGSSegment(
            ideal(R_oscar, e_polys),
            ideal(R_oscar, n_polys),
            ideal(R_oscar, b_polys),
        ))
    end
    return result
end

# Convert a Singular sideal to a Vector of Oscar polynomials.
function _sideal_to_oscar_vec(sing_ideal, SF, SR, R_oscar, params, main_vars)
    polys = Oscar.MPolyRingElem[]
    for i in 1:Oscar.Singular.ngens(sing_ideal)
        f_sing = sing_ideal[i]
        f_oscar = singular_to_oscar_poly(f_sing, SF, SR, R_oscar, params, main_vars)
        push!(polys, f_oscar)
    end
    return polys
end

# Parse P-representation segment description.
# seg_raw is a Vector{Any} where each element is [prime_sideal, [excl_sideals...]]
function _parse_segment(seg_raw, SF, SR, R_oscar, params, main_vars)
    segment = Vector[]
    for entry in seg_raw
        # entry[1] = prime sideal, entry[2] = Vector of excl sideals
        prime_sing = entry[1]
        excl_list_raw = entry[2]

        prime_polys = _sideal_to_oscar_vec(prime_sing, SF, SR, R_oscar, params, main_vars)
        prime_ideal = ideal(R_oscar, prime_polys)

        excl_ideals = Oscar.MPolyIdeal[]
        for excl_sing in excl_list_raw
            excl_polys = _sideal_to_oscar_vec(excl_sing, SF, SR, R_oscar, params, main_vars)
            push!(excl_ideals, ideal(R_oscar, excl_polys))
        end

        push!(segment, [prime_ideal, excl_ideals])
    end
    return segment
end
