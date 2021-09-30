##
using WavePropBase
using WavePropBase.Geometry
using WavePropBase.Trees
using StaticArrays
using Test

# recursively check that all point in a cluster tree are in the bounding box
function test_cluster_tree(clt)
    bbox = clt.container
    for iloc in clt.index_range
        x     = clt._elements[iloc]
        x ∈ bbox || (return false)
    end
    if !isroot(clt)
        clt ∈ clt.parent.children || (return false)
    end
    if !isleaf(clt)
        for child in clt.children
            test_cluster_tree(child) || (return false)
        end
    end
    return true
end

@testset "ClusterTree" begin
    @testset "1d" begin
        points    = SVector.([4,3,1,2,5,-1.0])
        splitter  = GeometricSplitter(nmax=1)
        clt = ClusterTree(points,splitter)
        @test sortperm(points) == clt.loc2glob
        splitter  = GeometricMinimalSplitter(nmax=1)
        clt = ClusterTree(points,splitter)
        @test sortperm(points) == clt.loc2glob
        splitter  = CardinalitySplitter(nmax=1)
        clt = ClusterTree(points,splitter)
        @test sortperm(points) == clt.loc2glob
        splitter  = PrincipalComponentSplitter(nmax=1)
        clt = ClusterTree(points,splitter)
        @test sortperm(points) == clt.loc2glob
        splitter  = DyadicSplitter(nmax=1)
        clt = ClusterTree(points,splitter)
        @test sortperm(points) == clt.loc2glob
    end

    @testset "2d" begin
        points    = rand(SVector{2,Float64},1000)
        splitter  = GeometricSplitter(nmax=1)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
        splitter  = GeometricMinimalSplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
        splitter   = CardinalitySplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
        splitter   = PrincipalComponentSplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
        splitter   = DyadicSplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
    end

    @testset "3d" begin
        points = rand(SVector{3,Float64},1000)
        splitter  = GeometricSplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
        splitter   = GeometricMinimalSplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
        splitter  = CardinalitySplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
        splitter   = PrincipalComponentSplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
        splitter   = DyadicSplitter(nmax=32)
        clt = ClusterTree(points,splitter)
        @test test_cluster_tree(clt)
    end
end
