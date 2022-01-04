using Test
using SafeTestsets

@safetestset "Polynomials" begin include("polynomials_test.jl") end

@safetestset "HyperRectangle" begin include("hyperrectangle_test.jl") end

@safetestset "Tensor lagrange interpolant" begin include("tensorlaginterp_test.jl") end

@safetestset "Lagrange Triangle" begin include("lagrangetriangle_test.jl") end
