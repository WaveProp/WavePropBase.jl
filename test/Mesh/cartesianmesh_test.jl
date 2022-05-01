using Test
import WavePropBase as WPB
using StaticArrays

@testset "One dimension" begin
    l = WPB.HyperRectangle(0.0,1.0)
    E = typeof(l) # type of mesh element
    h = 0.1 # grid spacing
    mesh = WPB.UniformCartesianMesh(l,(10,))
    @test mesh == WPB.UniformCartesianMesh(l,10)
    @test mesh == WPB.UniformCartesianMesh(WPB.grids(mesh)) == WPB.UniformCartesianMesh(WPB.grids(mesh)...)
    @test mesh == WPB.UniformCartesianMesh(l; step=step(mesh))
    @test WPB.xgrid(mesh) == LinRange(0.,1.,11)
    @test keys(mesh) == (E,)
    iter = WPB.ElementIterator(mesh,E)
    @test eltype(iter) == E
    @test iter[1] ≈ WPB.HyperRectangle(0,0.1)
    @test iter[3] ≈ WPB.HyperRectangle(0.2,0.3)
    for (n,el) in enumerate(iter)
        @test el ≈ WPB.HyperRectangle((n-1)*0.1,n*0.1)
    end
    iter = WPB.NodeIterator(mesh)
    @test eltype(iter) == SVector{1,Float64}
    @test size(iter) == (11,)
    for (n,x) in enumerate(iter)
        @test x ≈ SVector((n-1)*h)
    end
end

@testset "Two dimensions" begin
    l    = WPB.HyperRectangle((0.0,0.0),(1.0,1.0))
    E    = typeof(l) # type of mesh element
    mesh = WPB.UniformCartesianMesh(l,(10,20))
    iter = WPB.ElementIterator(mesh,E)
    @test mesh == WPB.UniformCartesianMesh(WPB.grids(mesh)) == WPB.UniformCartesianMesh(WPB.grids(mesh)...)
    @test mesh == WPB.UniformCartesianMesh(l; step=step(mesh))
    @test WPB.xgrid(mesh) == LinRange(0.,1.,11)
    @test WPB.ygrid(mesh) == LinRange(0.,1.,21)
    @test iter[1,1] ≈ WPB.HyperRectangle((0,0),(0.1,0.05))
    @test length(iter) == 200
    @test size(iter) == (10,20)
    for I in CartesianIndices(iter)
        i,j = Tuple(I)
        el = iter[i,j]
        @test el ≈ WPB.HyperRectangle(((i-1)*0.1,(j-1)*0.05),(i*0.1,j*0.05))
    end

    iter = WPB.NodeIterator(mesh)
    @test eltype(iter) == SVector{2,Float64}
    @test length(iter) == 11*21
    @test size(iter) == (11,21)
    for I in CartesianIndices(iter)
        i,j = Tuple(I)
        el = iter[i,j]
        @test el ≈ SVector((i-1)*0.1,(j-1)*0.05)
    end
end
