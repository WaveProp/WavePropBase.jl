"""
    module Mesh

Mesh data structures used to solve PDEs.
"""

module Mesh

using StaticArrays
using LinearAlgebra
using OrderedCollections

using WavePropBase.Utils
using WavePropBase.Geometry
using WavePropBase.Interpolation
using WavePropBase.Integration

# import all methods in WavePropBase.INTERFACE_METHODS
using WavePropBase
WavePropBase.@import_interface

export
    # abstract types
    AbstractMesh,
    # structs
    GenericMesh,
    ElementIterator,
    NodeIterator,
    SubMesh,
    UniformCartesianMesh,
    ParametricElement,
    # methods
    elements,
    ent2tags,
    el2qnodes,
    elements,
    near_interaction_list,
    vals,
    compute_quadrature!,
    dof,
    derivative,
    derivative2,
    measure,
    convert_to_2d,
    decompose,
    mesh,
    grids

include("abstractmesh.jl")
include("genericmesh.jl")
include("cartesianmesh.jl")
include("submesh.jl")
include("decompose.jl")

end
