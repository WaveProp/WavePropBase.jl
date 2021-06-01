using Test
using WavePropBase

@testset "One dimension" begin
    l = HyperRectangle(0.0,1.0)
    E = typeof(l) # type of mesh element
    h = 0.1
    mesh = CartesianMesh(domain=l,sz=(10,))
    keys(mesh) == (E,)
    iter = ElementIterator(mesh,E)
    @test iter[1] ≈ HyperRectangle(0,0.1)
    @test iter[3] ≈ HyperRectangle(0.2,0.3)
    for (n,el) in enumerate(iter)
        @test el ≈ HyperRectangle((n-1)*0.1,n*0.1)
    end
end

@testset "Two dimensions" begin
    l    = HyperRectangle((0.0,0.0),(1.0,1.0))
    E    = typeof(l) # type of mesh element
    mesh = CartesianMesh(domain=l,sz=(10,20))
    iter = ElementIterator(mesh,E)
    @test iter[1,1] ≈ HyperRectangle((0,0),(0.1,0.05))
    @test length(iter) == 200
    @test size(iter) == (10,20)
end
