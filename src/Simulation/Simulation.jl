module Simulation

using StaticArrays
using Printf
using RecipesBase
using OrderedCollections

using WavePropBase
using WavePropBase.Utils
using WavePropBase.Geometry
using WavePropBase.Integration
using WavePropBase.Interpolation
using WavePropBase.Mesh

WavePropBase.@import_interface

export
    # types
    AbstractPDE,
    Laplace,
    Helmholtz,
    Elastostatic,
    Maxwell

include("pde.jl")

end # module
