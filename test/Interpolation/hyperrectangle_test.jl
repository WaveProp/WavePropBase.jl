using Test
using WavePropBase
using WavePropBase.Geometry
using WavePropBase.Interpolation
using StaticArrays

@testset "HyperRectangle tests" begin
    low_corner  = SVector(0.0,0.0)
    high_corner = SVector(1.0,2.0)
    mid         = (low_corner + high_corner)/2
    rec = HyperRectangle(low_corner,high_corner)
    @test mid == center(rec)
    @test high_corner ∈ rec
    @test low_corner ∈ rec
    @test mid ∈ rec
    @test !in(high_corner + SVector(1,1),rec)
    rec1, rec2 = split(rec)
    @test low_corner ∈ rec1
    @test high_corner ∈ rec2
    @test !(low_corner ∈ rec2)
    @test !(high_corner ∈ rec1)
    @test diameter(rec) == sqrt(1^2 + 2^2)
    @test radius(rec) == sqrt(1^2 + 2^2)/2
    # bbox
    pts = SVector{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:1
            push!(pts,SVector(x,y))
        end
    end
    @test HyperRectangle(pts)      == HyperRectangle(SVector(-1.,-1),SVector(1,1.))
    @test HyperRectangle(pts,true) == HyperRectangle(SVector(-1.,-1),SVector(1,1.))
    pts = SVector{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:2
            push!(pts,SVector(x,y))
        end
    end
    @test HyperRectangle(pts)      == HyperRectangle(SVector(-1.,-1),SVector(1,2.))
    @test HyperRectangle(pts,true) == HyperRectangle(SVector(-1.5,-1),SVector(3/2,2.))
    rec1 = HyperRectangle(SVector(0,0),SVector(1,1))
    rec2 = HyperRectangle(SVector(2,0),SVector(3,1))
    @test distance(rec1,rec2) ≈ 1
    x  = SVector(0.5,0.5)
    @test distance(x,rec1) ≈ distance(rec1,x) ≈ 0
    x  = SVector(0.0,2.0)
    @test distance(x,rec1) ≈ distance(rec1,x) ≈ 1
    x  = SVector(2.0,2.0)
    @test distance(x,rec1) ≈ distance(rec1,x) ≈ √2
    rec2 = HyperRectangle(SVector(2,2),SVector(3,3))
    distance(rec1,rec2) ≈ sqrt(2)
end

@testset "HyperCube tests" begin
    low_corner  = SVector(1.5,-2.3)
    side = 2.2
    high_corner = low_corner .+ side
    mid         = (low_corner + high_corner)/2
    rec = HyperCube(low_corner,side)
    @test mid == center(rec)
    @test high_corner ∈ rec
    @test low_corner ∈ rec
    @test mid ∈ rec
    @test !in(high_corner + SVector(1,1),rec)
    rec1, rec2 = split(rec)
    @test low_corner ∈ rec1
    @test high_corner ∈ rec2
    @test !(low_corner ∈ rec2)
    @test !(high_corner ∈ rec1)
    @test diameter(rec) == side*sqrt(length(low_corner))
    @test radius(rec) == side*sqrt(length(low_corner))/2
    #=
    # bbox
    pts = SVector{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:1
            push!(pts,SVector(x,y))
        end
    end
    @test HyperRectangle(pts)      == HyperRectangle(SVector(-1.,-1),SVector(1,1.))
    @test HyperRectangle(pts,true) == HyperRectangle(SVector(-1.,-1),SVector(1,1.))
    pts = SVector{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:2
            push!(pts,SVector(x,y))
        end
    end
    @test HyperRectangle(pts)      == HyperRectangle(SVector(-1.,-1),SVector(1,2.))
    @test HyperRectangle(pts,true) == HyperRectangle(SVector(-1.5,-1),SVector(3/2,2.))
    =#
    rec1 = HyperCube(SVector(0,0),2)
    rec2 = HyperCube(SVector(3,0),1)
    @test distance(rec1,rec2) ≈ 1
    x  = SVector(0.5,0.5)
    @test distance(x,rec1) ≈ distance(rec1,x) ≈ 0
    x  = SVector(0.0,3.0)
    @test distance(x,rec1) ≈ distance(rec1,x) ≈ 1
    x  = SVector(3.0,3.0)
    @test distance(x,rec1) ≈ distance(rec1,x) ≈ √2
    rec2 = HyperCube(SVector(3,3),1)
    distance(rec1,rec2) ≈ sqrt(2)
end
