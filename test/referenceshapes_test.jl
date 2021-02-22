using Test
using WavePropBase

@testset "Line" begin
    l = ReferenceLine()
    @test ambient_dimension(l) == 1
    @test geometric_dimension(l) == 1
    x = 0.5
    @test x ∈ l
    x = 1.0
    @test x ∈ l
    x = 1.1
    @test !in(x,l)
end
@testset "Triangle" begin
    t = ReferenceTriangle()
    @test ambient_dimension(t) == 2
    @test geometric_dimension(t) == 2
    x = (0.5,0.5)
    @test x ∈ t
    x = (1.0,0.0)
    @test x ∈ t
    x = (1.1,0.0)
    @test !in(x,t)
end
@testset "Square" begin
    t = ReferenceSquare()
    @test ambient_dimension(t) == 2
    @test geometric_dimension(t) == 2
    x = (0.5,0.5)
    @test x ∈ t
    x = (1.0,0.0)
    @test x ∈ t
    x = (1.1,0.0)
    @test !in(x,t)
end
@testset "Tetrahedron" begin
    t = ReferenceTetrahedron()
    @test ambient_dimension(t) == 3
    @test geometric_dimension(t) == 3
    x = (0.5,0.5,0.0)
    @test x ∈ t
    x = (1.0,0.0,0.0) # point on edge
    @test x ∈ t
    x = (1.1,0.0,0.0)
    @test !in(x,t)
end
