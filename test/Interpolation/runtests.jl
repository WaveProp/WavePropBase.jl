using Test
using SafeTestsets

@safetestset "Polynomials" begin
    include("polynomials_test.jl")
end

@safetestset "HyperRectangle" begin
    include("hyperrectangle_test.jl")
end

@safetestset "Tensor lagrange interpolant" begin
    include("tensorlaginterp_test.jl")
end

@safetestset "Elements" begin
    include("element_test.jl")
end
