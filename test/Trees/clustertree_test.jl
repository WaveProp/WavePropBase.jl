using Test
import WavePropBase as WPB
using StaticArrays

# recursively check that all point in a cluster tree are in the bounding box
function test_cluster_tree(clt)
    bbox = WPB.container(clt)
    for iloc in WPB.index_range(clt)
        x = WPB.root_elements(clt)[iloc]
        x ∈ bbox || (return false)
    end
    if !WPB.isroot(clt)
        clt ∈ clt.parent.children || (return false)
    end
    if !WPB.isleaf(clt)
        for child in clt.children
            test_cluster_tree(child) || (return false)
        end
    end
    return true
end

@testset "ClusterTree" begin
    @testset "1d" begin
        points = SVector.([4, 3, 1, 2, 5, -1.0])
        splitter = WPB.GeometricSplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test sortperm(points) == clt.loc2glob
        splitter = WPB.GeometricMinimalSplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test sortperm(points) == clt.loc2glob
        splitter = WPB.CardinalitySplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test sortperm(points) == clt.loc2glob
        splitter = WPB.PrincipalComponentSplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test sortperm(points) == clt.loc2glob
        splitter = WPB.DyadicSplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test sortperm(points) == clt.loc2glob
        splitter = WPB.DyadicMinimalSplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test sortperm(points) == clt.loc2glob
    end

    @testset "2d" begin
        points = rand(SVector{2,Float64}, 1000)
        splitter = WPB.GeometricSplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.GeometricMinimalSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.CardinalitySplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.PrincipalComponentSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.DyadicSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.DyadicMinimalSplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
    end

    @testset "3d" begin
        points = rand(SVector{3,Float64}, 1000)
        splitter = WPB.GeometricSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.GeometricMinimalSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.CardinalitySplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.PrincipalComponentSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.DyadicSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
        splitter = WPB.DyadicMinimalSplitter(; nmax=1)
        clt = WPB.ClusterTree(points, splitter)
        @test test_cluster_tree(clt)
    end

    @testset "3d + threads" begin
        threads = true
        points = rand(SVector{3,Float64}, 1000)
        splitter = WPB.GeometricSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter; threads)
        @test test_cluster_tree(clt)
        splitter = WPB.GeometricMinimalSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter; threads)
        @test test_cluster_tree(clt)
        splitter = WPB.CardinalitySplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter; threads)
        @test test_cluster_tree(clt)
        splitter = WPB.PrincipalComponentSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter; threads)
        @test test_cluster_tree(clt)
        splitter = WPB.DyadicSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter; threads)
        @test test_cluster_tree(clt)
        splitter = WPB.DyadicMinimalSplitter(; nmax=32)
        clt = WPB.ClusterTree(points, splitter; threads)
        @test test_cluster_tree(clt)
    end
end
