using Test
using WavePropBase
using StaticArrays

# TODO: complete these tests
@testset "Lagrange elements" begin
    @testset "LagrangeElement{ReferenceLine}" begin
        ## line in 3d
        vtx = ((0.,0.,0.),(1.,1.,1.))
        l   = LagrangeElement{ReferenceLine}(vtx)
        @test domain(l) == ReferenceLine()
        ## line in 2d
        a = (0.,0.)
        b = (1.,1.)
        l   = LagrangeElement{ReferenceLine}((a,b))
        @test domain(l) == ReferenceLine()
    end
    @testset "LagrangeElement{ReferenceTriangle}" begin
        # triangle in 2d
        vtx = SVector(SVector(0.,0.),SVector(0.,1.),SVector(-1.,0))
        t   = LagrangeElement{ReferenceTriangle}(vtx)
        # triangle in 3d
        vtx = SVector(SVector(0.,0.,0.),SVector(0.,1.,0.),SVector(-1.,0,0.))
        t   = LagrangeElement{ReferenceTriangle}(vtx)
    end
    @testset "Tetrahedron" begin
    end
end
