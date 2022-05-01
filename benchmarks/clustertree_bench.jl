import WavePropBase as WPB
using StaticArrays
using BenchmarkTools

n = 1_000_000

pts      = rand(WPB.Point3D,n)
splitter = WPB.GeometricSplitter(nmax=100)
@btime WPB.ClusterTree($pts,$splitter;threads=false)
