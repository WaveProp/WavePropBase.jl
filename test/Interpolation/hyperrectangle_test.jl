using Test
import WavePropBase as WPB
using StaticArrays

@testset "HyperRectangle tests" begin
    low_corner  = SVector(0.0,0.0)
    high_corner = SVector(1.0,2.0)
    mid         = (low_corner + high_corner)/2
    rec = WPB.HyperRectangle(low_corner,high_corner)
    @test mid == WPB.center(rec)
    @test high_corner ∈ rec
    @test low_corner ∈ rec
    @test mid ∈ rec
    @test !in(high_corner + SVector(1,1),rec)
    rec1, rec2 = split(rec)
    @test low_corner ∈ rec1
    @test high_corner ∈ rec2
    @test !(low_corner ∈ rec2)
    @test !(high_corner ∈ rec1)
    @test WPB.diameter(rec) == sqrt(1^2 + 2^2)
    @test WPB.radius(rec) == sqrt(1^2 + 2^2)/2
    # bbox
    pts = SVector{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:1
            push!(pts,SVector(x,y))
        end
    end
    @test WPB.HyperRectangle(pts)      == WPB.HyperRectangle(SVector(-1.,-1),SVector(1,1.))
    @test WPB.HyperRectangle(pts,true) == WPB.HyperRectangle(SVector(-1.,-1),SVector(1,1.))
    pts = SVector{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:2
            push!(pts,SVector(x,y))
        end
    end
    @test WPB.HyperRectangle(pts)      == WPB.HyperRectangle(SVector(-1.,-1),SVector(1,2.))
    @test WPB.HyperRectangle(pts,true) == WPB.HyperRectangle(SVector(-1.5,-1),SVector(3/2,2.))
    rec1 = WPB.HyperRectangle(SVector(0,0),SVector(1,1))
    rec2 = WPB.HyperRectangle(SVector(2,0),SVector(3,1))
    @test WPB.distance(rec1,rec2) ≈ 1
    x  = SVector(0.5,0.5)
    @test WPB.distance(x,rec1) ≈ WPB.distance(rec1,x) ≈ 0
    x  = SVector(0.0,2.0)
    @test WPB.distance(x,rec1) ≈ WPB.distance(rec1,x) ≈ 1
    x  = SVector(2.0,2.0)
    @test WPB.distance(x,rec1) ≈ WPB.distance(rec1,x) ≈ √2
    rec2 = WPB.HyperRectangle(SVector(2,2),SVector(3,3))
    WPB.distance(rec1,rec2) ≈ sqrt(2)
end

@testset "HyperCube tests" begin
    low_corner  = SVector(1.5,-2.3)
    side = 2.2
    high_corner = low_corner .+ side
    mid         = (low_corner + high_corner)/2
    rec = WPB.HyperCube(low_corner,side)
    @test mid == WPB.center(rec)
    @test high_corner ∈ rec
    @test low_corner ∈ rec
    @test mid ∈ rec
    @test !in(high_corner + SVector(1,1),rec)
    rec1, rec2 = split(rec)
    @test low_corner ∈ rec1
    @test high_corner ∈ rec2
    @test !(low_corner ∈ rec2)
    @test !(high_corner ∈ rec1)
    @test WPB.diameter(rec) == side*sqrt(length(low_corner))
    @test WPB.radius(rec) == side*sqrt(length(low_corner))/2
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
    rec1 = WPB.HyperCube(SVector(0,0),2)
    rec2 = WPB.HyperCube(SVector(3,0),1)
    @test WPB.distance(rec1,rec2) ≈ 1
    x  = SVector(0.5,0.5)
    @test WPB.distance(x,rec1) ≈ WPB.distance(rec1,x) ≈ 0
    x  = SVector(0.0,3.0)
    @test WPB.distance(x,rec1) ≈ WPB.distance(rec1,x) ≈ 1
    x  = SVector(3.0,3.0)
    @test WPB.distance(x,rec1) ≈ WPB.distance(rec1,x) ≈ √2
    rec2 = WPB.HyperCube(SVector(3,3),1)
    WPB.distance(rec1,rec2) ≈ sqrt(2)
end
