using SafeTestsets

@safetestset "Reference shapes" begin include("referenceshapes_test.jl") end

@safetestset "HyperRectangle" begin include("hyperrectangle_test.jl") end

@safetestset "Domain" begin include("domain_test.jl") end

@safetestset "Lagrange interpolation" begin include("lagrangeinterp_test.jl") end

@safetestset "Elements" begin include("element_test.jl") end

@safetestset "Cartesian mesh" begin include("cartesianmesh_test.jl") end
