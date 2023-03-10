using Test
import WavePropBase as WPB
using StaticArrays

@test WPB.svector(i->i^2,3) == SVector(ntuple(i->i^2, 3))

T = SMatrix{2,2,Float64,4}

A = WPB.PseudoBlockMatrix{T}(undef, 2, 2)
