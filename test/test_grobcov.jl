@testset "grobcov" verbose = true begin

    # -----------------------------------------------------------------------
    # Single-parameter examples
    # -----------------------------------------------------------------------

    @testset "I = [a*x - 1]: 1 param, 1 variable" begin
        R, (a, x) = polynomial_ring(QQ, ["a", "x"])
        result = grobcov(ideal(R, [a * x - 1]), [a], [x])

        # Groebner cover has exactly 2 segments:
        #   generic (a ≠ 0): lpp = [x],      basis = {a*x - 1}
        #   special (a = 0): lpp = [1],      basis = {1}
        @test result isa Vector{GrobCovSegment}
        @test length(result) == 2

        generic_seg, special_seg = result

        # Generic segment
        @test generic_seg.lpp == [x]
        @test a * x - 1 ∈ gens(generic_seg.basis)
        @test all(p -> parent(p) === R, generic_seg.lpp)
        @test all(p -> parent(p) === R, gens(generic_seg.basis))

        # Special segment (degenerate case)
        @test special_seg.lpp == [one(R)]
        @test one(R) ∈ gens(special_seg.basis)
    end

    @testset "I = [a*x - 1, a*y - 1]: 1 param, 2 variables" begin
        R, (a, x, y) = polynomial_ring(QQ, ["a", "x", "y"])
        result = grobcov(ideal(R, [a * x - 1, a * y - 1]), [a], [x, y])

        @test length(result) == 2
        generic_seg = result[1]

        # Both x and y appear as leading power products on the generic segment
        @test issubset([x, y], generic_seg.lpp)
        @test all(p -> parent(p) === R, gens(generic_seg.basis))
    end

    # -----------------------------------------------------------------------
    # Two-parameter example
    # -----------------------------------------------------------------------

    @testset "I = [a*x^2 + b*y, x - b*y]: 2 params, 2 variables" begin
        R, (a, b, x, y) = polynomial_ring(QQ, ["a", "b", "x", "y"])
        result = grobcov(ideal(R, [a * x^2 + b * y, x - b * y]), [a, b], [x, y])

        @test result isa Vector{GrobCovSegment}
        @test length(result) == 3

        for seg in result
            @test seg.lpp isa Vector{<:Oscar.MPolyRingElem}
            @test all(p -> parent(p) === R, seg.lpp)
            @test all(p -> parent(p) === R, gens(seg.basis))
        end
    end

    # -----------------------------------------------------------------------
    # P-representation structure
    # -----------------------------------------------------------------------

    @testset "segment P-representation is well-formed" begin
        R, (a, x) = polynomial_ring(QQ, ["a", "x"])
        result = grobcov(ideal(R, [a * x - 1]), [a], [x])

        for seg in result
            @test seg.segment isa Vector
            for entry in seg.segment
                @test entry isa Vector
                @test length(entry) == 2
                prime, excls = entry
                @test prime isa Oscar.MPolyIdeal
                @test excls isa Vector{<:Oscar.MPolyIdeal}
                @test all(p -> parent(p) === R, gens(prime))
                for excl in excls
                    @test all(p -> parent(p) === R, gens(excl))
                end
            end
        end

        # Generic segment: V(0) \ V(a)  →  prime = 0,  excl = [ideal(a)]
        prime, excls = result[1].segment[1]
        @test iszero(only(gens(prime)))
        @test a ∈ gens(only(excls))

        # Special segment: V(a) \ V(1)  →  prime = a,  excl = [ideal(1)]
        prime2, excls2 = result[2].segment[1]
        @test a ∈ gens(prime2)
    end

end

@testset "cgsdr" verbose = true begin

    # -----------------------------------------------------------------------
    # Single-parameter examples
    # -----------------------------------------------------------------------

    @testset "I = [a*x - 1]: 1 param, 1 variable" begin
        R, (a, x) = polynomial_ring(QQ, ["a", "x"])
        result = cgsdr(ideal(R, [a * x - 1]), [a], [x])

        # 2 segments:
        #   seg 1: null = 0  (full space),  nonnull = a  (a ≠ 0),  basis = {a*x - 1}
        #   seg 2: null = a  (a = 0),       nonnull = 1,            basis = {1}
        @test result isa Vector{CGSSegment}
        @test length(result) == 2

        seg1, seg2 = result

        @test iszero(only(gens(seg1.null)))
        @test a     ∈ gens(seg1.nonnull)
        @test a * x - 1 ∈ gens(seg1.basis)

        @test a      ∈ gens(seg2.null)
        @test one(R) ∈ gens(seg2.basis)

        @test all(p -> parent(p) === R, gens(seg1.null))
        @test all(p -> parent(p) === R, gens(seg1.nonnull))
        @test all(p -> parent(p) === R, gens(seg1.basis))
        @test all(p -> parent(p) === R, gens(seg2.null))
        @test all(p -> parent(p) === R, gens(seg2.basis))
    end

    @testset "I = [a*x - 1, a*y - 1]: 1 param, 2 variables" begin
        R, (a, x, y) = polynomial_ring(QQ, ["a", "x", "y"])
        result = cgsdr(ideal(R, [a * x - 1, a * y - 1]), [a], [x, y])

        @test result isa Vector{CGSSegment}
        @test length(result) >= 1
        for seg in result
            @test all(p -> parent(p) === R, gens(seg.null))
            @test all(p -> parent(p) === R, gens(seg.nonnull))
            @test all(p -> parent(p) === R, gens(seg.basis))
        end
    end

    @testset "each segment field is MPolyIdeal" begin
        R, (a, b, x, y) = polynomial_ring(QQ, ["a", "b", "x", "y"])
        result = cgsdr(ideal(R, [a * x^2 + b * y, x - b * y]), [a, b], [x, y])

        for seg in result
            @test seg.null    isa Oscar.MPolyIdeal
            @test seg.nonnull isa Oscar.MPolyIdeal
            @test seg.basis   isa Oscar.MPolyIdeal
        end
    end

end
