"""
    module Trees

Module defining some commonly used tree structures, as well as an interface for
talking about trees inside the `WaveProp` organization.
"""
module Trees

import Base: ==, length, in, iterate, getindex, lastindex, isempty, eltype, keys
import Base: union, setdiff, intersect, issubset

using StaticArrays
using LinearAlgebra
using AbstractTrees
using Printf
using Statistics: median

using WavePropBase
using WavePropBase.Utils
using WavePropBase.Geometry
using WavePropBase.Interpolation

WavePropBase.@import_interface

include("abstracttree.jl")
include("clustertree.jl")
include("splitter.jl")

export
    # abstract types
    AbstractTree,
    # types
    ClusterTree,
    GeometricSplitter,
    GeometricMinimalSplitter,
    PrincipalComponentSplitter,
    DyadicSplitter,
    CardinalitySplitter,
    # re-exported from AbstractTrees
    Leaves,
    PreOrderDFS,
    PostOrderDFS,
    print_tree,
    # functions
    children,
    parent,
    isroot,
    isleaf,
    depth,
    container,
    loc2glob,
    root_elements,
    index_range,
    filter_tree,
    filter_tree!

end # module
