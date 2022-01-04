using Test
using WavePropBase.Geometry
using WavePropBase.Interpolation
using LinearAlgebra
using ForwardDiff
using StaticArrays
import Random
Random.seed!(1)

# TODO: combine this file with element_test.jl
@testset "LagrangeTriangle" begin
    # triangle in 3d
    vtx = SVector(SVector(0.,0.,0.),SVector(0.,1.,0.),SVector(-1.,0,0.))
    t   = LagrangeTriangle(vtx)
    @test Interpolation.normal(t,SVector(0.1,0.1)) == SVector(0,0,1.)
    # lagrange nodes in parametric space, following Gmsh convention
    T = SVector{2,Float64}
    triangle3_nodes = [T(0,0),T(1,0),T(0,1)]
    triangle6_nodes = [T(0,0),T(1,0),T(0,1),T(1/2,0),T(1/2,1/2),T(0,1/2)]
    triangle10_nodes = [T(0,0),T(1,0),T(0,1),T(1/3,0),T(2/3,0),T(2/3,1/3),
                        T(1/3,2/3),T(0,2/3),T(0,1/3),T(1/3,1/3)]
    triangle_list = [triangle3_nodes,triangle6_nodes,triangle10_nodes]
    p = T(2/5,3/7)   # a point in the Reference Triangle
    @test p ∈ ReferenceTriangle()
    for triangle_nodes in triangle_list
        Np = length(triangle_nodes)
        # generate points in 3D space
        vtx = @SVector rand(SVector{3,Float64},Np)
        t = LagrangeTriangle(vtx)
        # test lagrange basis at lagrange nodes
        for (u,v) in zip(triangle_nodes,vtx)
            @test t(u) ≈ v
        end
        # test derivative
        @test Interpolation.jacobian(t,p) ≈ ForwardDiff.jacobian(t,p)
    end
end