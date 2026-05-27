using Test
using Oscar
using GrobCovBridge

@testset "GrobCovBridge" verbose = true begin
    include("test_conversion.jl")
    include("test_grobcov.jl")
end
