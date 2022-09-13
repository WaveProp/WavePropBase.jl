using SafeTestsets

@safetestset "Geometry" begin include("Geometry/runtests.jl") end

@safetestset "Interpolation" begin include("Interpolation/runtests.jl") end

@safetestset "Integration" begin include("Integration/runtests.jl") end

@safetestset "Mesh" begin include("Mesh/runtests.jl") end

@safetestset "ParametricEntities" begin include("ParametricEntities/runtests.jl") end

@safetestset "Trees" begin include("Trees/runtests.jl") end

# FIXME: I still can't figure out how to add an unregistered package as a
# dependency of the tests. I get some "unable to merge project" error from Pkg.
# Until that is resolved, the tests below are commented for the CI not to fail.
# @safetestset "IO" begin include("IO/runtests.jl") end
