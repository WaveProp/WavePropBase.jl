"""
    module Geometry

Module defining basic geometrical concepts.
"""
module Geometry

import Base: ==, length, in, iterate, getindex, lastindex, isempty, eltype, keys
import Base: union, setdiff, intersect, issubset

using StaticArrays
using LinearAlgebra
using OrderedCollections
using AbstractTrees
using RecipesBase
using Printf
using Statistics: median

using WavePropBase
using WavePropBase.Utils

WavePropBase.@import_interface

include("hyperrectangle.jl")
include("clustertree.jl")
include("referenceshapes.jl")
include("entities.jl")
include("domain.jl")

export
    # abstract types
    AbstractReferenceShape,
    AbstractEntity,
    AbstractParametricBody,
    # types
    ElementaryEntity,
    ParametricEntity,
    ParametricCurve,
    Domain,
    ReferencePoint,
    ReferenceLine,
    ReferenceTriangle,
    ReferenceTetrahedron,
    ReferenceSquare,
    HyperCube,
    HyperRectangle,
    ClusterTree,
    GeometricSplitter,
    GeometricMinimalSplitter,
    PrincipalComponentSplitter,
    DyadicSplitter,
    CardinalitySplitter,
    # functions
    clear_entities!,
    tag,
    key,
    assertequaldim,
    boundary,
    skeleton,
    internal_boundary,
    external_boundary,
    low_corner,
    high_corner,
    measure,
    vertices,
    line,
    new_tag,
    global_add_entity!,
    # global variables
    TAGS,
    ENTITIES

end # module
