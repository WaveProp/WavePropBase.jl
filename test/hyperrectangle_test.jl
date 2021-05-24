using Test
using WavePropBase
using StaticArrays

const WPBase = WavePropBase

@testset "Basic tests" begin
    low_corner  = (0.0,0.0)
    high_corner = (1.0,2.0)
    mid         = (low_corner .+ high_corner) ./ 2 |> SVector
    rec = WPBase.HyperRectangle(low_corner,high_corner)
    @test mid == WPBase.center(rec)
    @test (mid ∈ rec)  == true
    @test !in(high_corner .+ (1,1),rec)
    @test WPBase.diameter(rec) == sqrt(1^2 + 2^2)
    @test WPBase.radius(rec) == sqrt(1^2 + 2^2)/2
    # bbox
    pts = SVector{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:1
            push!(pts,(x,y))
        end
    end
    @test WPBase.bounding_box(pts) == WPBase.HyperRectangle((-1.,-1),(1,1.))
end
