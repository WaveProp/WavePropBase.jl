using Test
using SafeTestsets

@safetestset "Quadrature rules" begin include("quadrule_test.jl") end
