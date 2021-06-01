using Test
using WavePropBase

@testset "One dimension" begin
    l = HyperRectangle(0.0,1.0)
    h = 0.1
    mesh = CartesianMesh(domain=l,sz=(10,))
    @test mesh[1] ≈ HyperRectangle(0,0.1)
    @test mesh[3] ≈ HyperRectangle(0.2,0.3)
    for (n,el) in enumerate(mesh)
        @test el ≈ HyperRectangle((n-1)*0.1,n*0.1)
    end
    E = typeof(l) # type of mesh element
    keys(mesh) == (E,)
    iter = ElementIterator(mesh,E)
    for (n,el) in enumerate(mesh)
        @test el ≈ HyperRectangle((n-1)*0.1,n*0.1)
    end
end

@testset "Two dimensions" begin
    l = HyperRectangle((0.0,0.0),(1.0,1.0))
    mesh = CartesianMesh(domain=l,sz=(10,20))
    @test mesh[1,1] ≈ HyperRectangle((0,0),(0.1,0.05))
    @test length(mesh) == 200
    @test size(mesh) == (10,20)
    E    = typeof(l) # type of mesh element
    iter = ElementIterator(mesh,E)
end
