import Oscar.Singular as S

# ---------------------------------------------------------------------------
# Oscar → Singular
# ---------------------------------------------------------------------------

"""
    oscar_to_singular_parametric(f, R, params, main_vars)

Convert an Oscar polynomial `f` in ring `R` to a Singular polynomial in a
parametric ring `Q(params)[main_vars]`.

Returns `(SF, SR, sing_f)`:
- `SF`:     Singular `N_FField` (the parameter function field)
- `SR`:     Singular `PolyRing` over `SF`
- `sing_f`: `f` as an element of `SR`
"""
function oscar_to_singular_parametric(f, R, params, main_vars)
    _check_var_partition(R, params, main_vars)
    param_names = string.(symbols(R)[_var_indices(R, params)])
    main_names  = string.(symbols(R)[_var_indices(R, main_vars)])

    SF, sa_vars = S.FunctionField(S.QQ, param_names)
    SR, sx_vars = S.polynomial_ring(SF, main_names)

    param_indices = _var_indices(R, params)
    main_indices  = _var_indices(R, main_vars)

    sing_f = _oscar_poly_to_sing(f, SF, SR, sa_vars, sx_vars, param_indices, main_indices)
    return SF, SR, sing_f
end

"""
    oscar_ideal_to_singular_parametric(I, params, main_vars)

Convert an Oscar `MPolyIdeal` to a Singular `sideal` in a parametric ring.

Returns `(SF, SR, sing_I)`.
"""
function oscar_ideal_to_singular_parametric(I, params, main_vars)
    R = base_ring(I)
    gs = gens(I)
    isempty(gs) && error("Ideal must have at least one generator.")

    _check_var_partition(R, params, main_vars)
    param_names = string.(symbols(R)[_var_indices(R, params)])
    main_names  = string.(symbols(R)[_var_indices(R, main_vars)])

    SF, sa_vars = S.FunctionField(S.QQ, param_names)
    SR, sx_vars = S.polynomial_ring(SF, main_names)

    param_indices = _var_indices(R, params)
    main_indices  = _var_indices(R, main_vars)

    sing_polys = [_oscar_poly_to_sing(f, SF, SR, sa_vars, sx_vars, param_indices, main_indices)
                  for f in gs]
    sing_I = S.Ideal(SR, sing_polys)
    return SF, SR, sing_I
end

# Convert a single Oscar polynomial to a Singular polynomial in SR.
function _oscar_poly_to_sing(f, SF, SR, sa_vars, sx_vars, param_indices, main_indices)
    result = zero(SR)
    for (coeff, mono) in zip(coefficients(f), monomials(f))
        exp_vec = exponent_vector(mono, 1)

        # Build the parameter part as an element of SF
        param_part = SF(QQ(coeff))
        for (pi, si) in zip(param_indices, sa_vars)
            e = exp_vec[pi]
            if e > 0
                param_part *= si^e
            end
        end

        # Build the main-variable monomial in SR
        main_mono = one(SR)
        for (mi, sv) in zip(main_indices, sx_vars)
            e = exp_vec[mi]
            if e > 0
                main_mono *= sv^e
            end
        end

        result += SR(param_part) * main_mono
    end
    return result
end

# Return the (1-based) indices of `vars` in ring `R`.
function _var_indices(R, vars)
    all_gens = gens(R)
    return [findfirst(==(v), all_gens) for v in vars]
end

# Validate that `params` and `main_vars` form a partition of the generators of `R`.
# Without this, any generator that is neither a parameter nor a main variable would
# be silently ignored during conversion, corrupting the polynomial.
function _check_var_partition(R, params, main_vars)
    all_gens = gens(R)
    declared = vcat(collect(params), collect(main_vars))

    for v in declared
        findfirst(==(v), all_gens) === nothing &&
            throw(ArgumentError("Variable $v is not a generator of the base ring."))
    end
    for p in params, m in main_vars
        p == m && throw(ArgumentError("Variable $p appears in both params and main_vars."))
    end
    for g in all_gens
        findfirst(==(g), declared) === nothing &&
            throw(ArgumentError(
                "Variable $g is neither a parameter nor a main variable; params and " *
                "main_vars must together cover all variables of the base ring."))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Singular → Oscar
# ---------------------------------------------------------------------------

"""
    singular_to_oscar_poly(f_sing, SF, SR, R_oscar, params, main_vars)

Convert a Singular polynomial `f_sing` (element of `SR`, a PolyRing over an
`N_FField SF`) back to an Oscar polynomial in `R_oscar`.

Coefficients of `f_sing` must be polynomials in the parameters (i.e. rational
functions with a constant denominator). A non-constant denominator — a genuine
parametric rational function — cannot be represented in `R_oscar` and raises an
error rather than silently dropping the denominator.
"""
function singular_to_oscar_poly(f_sing, SF, SR, R_oscar, params, main_vars)
    SF_param_names = [string(S.symbols(SF)[i]) for i in 1:S.transcendence_degree(SF)]
    SR_param, sp_vars = S.polynomial_ring(S.QQ, SF_param_names)

    result = zero(R_oscar)

    for (c, m) in zip(S.coefficients(f_sing), S.monomials(f_sing))
        # c :: n_transExt, m :: spoly (monomial in main vars)
        # Get exponent vector for main variables
        main_expvec = [S.degree(m, i) for i in 1:S.nvars(SR)]

        # Convert numerator and denominator of c to polynomials in SR_param
        c_num = S.n_transExt_to_spoly(S.numerator(c); parent_ring = SR_param)
        c_den = S.n_transExt_to_spoly(S.denominator(c); parent_ring = SR_param)

        if !S.is_constant(c_den)
            error("singular_to_oscar_poly: coefficient has a non-constant denominator " *
                  "(parametric rational function). Only polynomial coefficients in the " *
                  "parameters are supported. Coefficient: $c")
        end
        den_coeff = _sing_constant_to_QQ(c_den, SR_param)

        for (cc, cm) in zip(S.coefficients(c_num), S.monomials(c_num))
            # cc :: n_Q (rational number), cm :: monomial in SR_param
            rat = QQ(S.QQ(cc))
            if den_coeff !== nothing
                rat = rat // den_coeff
            end
            param_expvec = [S.degree(cm, i) for i in 1:S.nvars(SR_param)]

            # Build the full exponent vector in R_oscar: params first, then main vars
            param_indices = _var_indices(R_oscar, params)
            main_indices  = _var_indices(R_oscar, main_vars)

            term = R_oscar(rat)
            for (pi, pe) in zip(param_indices, param_expvec)
                if pe > 0
                    term *= gen(R_oscar, pi)^pe
                end
            end
            for (mi, me) in zip(main_indices, main_expvec)
                if me > 0
                    term *= gen(R_oscar, mi)^me
                end
            end
            result += term
        end
    end
    return result
end

# If sing_poly is a scalar constant, return its QQ value; otherwise return nothing.
function _sing_constant_to_QQ(p, SR_param)
    S.is_constant(p) || return nothing
    S.iszero(p) && return QQ(0)
    # p is a nonzero constant: it has exactly one term
    coeffs = collect(S.coefficients(p))
    isempty(coeffs) && return QQ(0)
    return QQ(S.QQ(first(coeffs)))
end
