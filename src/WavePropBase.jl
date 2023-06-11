module WavePropBase

const PROJECT_ROOT = pkgdir(WavePropBase)

using StaticArrays
using LinearAlgebra
using SparseArrays
using NearestNeighbors
using DataStructures
using AbstractTrees
using Requires
using RecipesBase
using Printf
using Statistics: median
using SpecialFunctions
using HCubature
using WriteVTK
using ForwardDiff
using Gmsh
using Lebedev

include("Utils/utils.jl")
include("Geometry/Geometry.jl")
include("Interpolation/Interpolation.jl")
include("Integration/Integration.jl")
include("Mesh/Mesh.jl")
include("IntegralEquations/IntegralEquations.jl")
include("Utils/recipes.jl")

function __init__()
    # handling of optional dependencies
    # @require WriteVTK = "64499a7a-5c06-52f2-abe2-ccb03c286192" include("Utils/writevtk_ext.jl")
    # @require Gmsh = "705231aa-382f-11e9-3f0c-b7cb4346fdeb" include("Utils/gmsh_ext.jl")
    # @require ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210" include("Utils/forwarddiff_ext.jl")
end

end
