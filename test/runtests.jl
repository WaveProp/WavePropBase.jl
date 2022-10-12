using SafeTestsets

@safetestset "Geometry" begin include("Geometry/runtests.jl") end

@safetestset "Interpolation" begin include("Interpolation/runtests.jl") end

@safetestset "Integration" begin include("Integration/runtests.jl") end

@safetestset "Mesh" begin include("Mesh/runtests.jl") end

@safetestset "ParametricEntities" begin include("ParametricEntities/runtests.jl") end

@safetestset "Trees" begin include("Trees/runtests.jl") end

@safetestset "IO" begin include("IO/runtests.jl") end
