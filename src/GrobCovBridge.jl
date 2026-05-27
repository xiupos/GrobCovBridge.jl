module GrobCovBridge

using Oscar

export GrobCovSegment, CGSSegment
export grobcov, cgsdr

include("conversion.jl")
include("output.jl")
include("grobcov.jl")

end
