@testset "Conversion" verbose = true begin

    # -----------------------------------------------------------------------
    # Oscar → Singular direction
    # -----------------------------------------------------------------------

    @testset "oscar_to_singular_parametric: ring types and dimensions" begin
        R, (a, b, x, y) = polynomial_ring(QQ, ["a", "b", "x", "y"])
        SF, SR, _ = GrobCovBridge.oscar_to_singular_parametric(a * x, R, [a, b], [x, y])

        @test SF isa Oscar.Singular.N_FField
        @test SR isa Oscar.Singular.PolyRing{<:Oscar.Singular.n_transExt}
        @test Oscar.Singular.transcendence_degree(SF) == 2
        @test Oscar.Singular.nvars(SR) == 2
    end

    @testset "oscar_to_singular_parametric: output is spoly" begin
        R, (a, b, x, y) = polynomial_ring(QQ, ["a", "b", "x", "y"])
        cases = [a * x, a^2 - b, x^2 + y^2 - 1, a^2 * x + 2 * a * b * y - b^3]
        @testset "f = $(repr(f))" for f in cases
            _, _, sing_f = GrobCovBridge.oscar_to_singular_parametric(f, R, [a, b], [x, y])
            @test sing_f isa Oscar.Singular.spoly
        end
    end

    @testset "oscar_ideal_to_singular_parametric: generator count preserved" begin
        R, (a, x, y) = polynomial_ring(QQ, ["a", "x", "y"])
        for n in [1, 2, 3]
            gs = [a^i * x^(n - i) - i for i in 0:(n - 1)]
            _, _, sing_I = GrobCovBridge.oscar_ideal_to_singular_parametric(
                ideal(R, gs), [a], [x, y])
            @test sing_I isa Oscar.Singular.sideal
            @test Oscar.Singular.ngens(sing_I) == n
        end
    end

    # -----------------------------------------------------------------------
    # Singular → Oscar direction
    # -----------------------------------------------------------------------

    @testset "singular_to_oscar_poly: concrete values" begin
        R, (a, x, y) = polynomial_ring(QQ, ["a", "x", "y"])
        SF, (sa,) = Oscar.Singular.FunctionField(Oscar.Singular.QQ, ["a"])
        SR, (sx, sy) = Oscar.Singular.polynomial_ring(SF, ["x", "y"])

        cases = [
            (SR(sa) * sx + one(SR),         a * x + 1),
            (SR(sa^2) - 2 * one(SR),         a^2 - 2),
            (SR(sa) * sy - 3 * sx,           a * y - 3 * x),
            (one(SR),                         one(R)),
        ]
        @testset "f_sing → $(repr(expected))" for (f_sing, expected) in cases
            result = GrobCovBridge.singular_to_oscar_poly(f_sing, SF, SR, R, [a], [x, y])
            @test result == expected
        end
    end

    # -----------------------------------------------------------------------
    # Round-trip tests
    # -----------------------------------------------------------------------

    @testset "round-trip Oscar → Singular → Oscar" begin
        R, (a, b, x, y) = polynomial_ring(QQ, ["a", "b", "x", "y"])
        params    = [a, b]
        main_vars = [x, y]

        cases = [
            a * x + b * y,
            a^2 * x - b^2 * y + 1,
            3 * a * b * x^2 + x - 1,
            a + b,          # no main variables
            R(QQ(3, 2)),    # rational constant
        ]
        @testset "f = $(repr(f))" for f in cases
            SF, SR, f_sing = GrobCovBridge.oscar_to_singular_parametric(
                f, R, params, main_vars)
            @test GrobCovBridge.singular_to_oscar_poly(
                f_sing, SF, SR, R, params, main_vars) == f
        end
    end

    # -----------------------------------------------------------------------
    # Error handling
    # -----------------------------------------------------------------------

    @testset "empty params list throws ArgumentError" begin
        R, (x,) = polynomial_ring(QQ, ["x"])
        @test_throws ArgumentError GrobCovBridge.oscar_ideal_to_singular_parametric(
            ideal(R, [x - 1]), QQMPolyRingElem[], [x])
    end

end
