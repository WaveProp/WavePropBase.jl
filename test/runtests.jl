using WavePropBase
using SafeTestsets

@safetestset "Reference shapes" begin include("referenceshapes_test.jl") end

@safetestset "HyperRectangle" begin include("hyperrectangle_test.jl") end
