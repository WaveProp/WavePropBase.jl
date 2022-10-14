module WavePropBase

const PROJECT_ROOT =  pkgdir(WavePropBase)

using StaticArrays
using LinearAlgebra
using AbstractTrees
using RecipesBase
using Requires
using Printf
using Statistics: median

include("Utils/Utils.jl")
include("Geometry/Geometry.jl")
include("Interpolation/Interpolation.jl")
include("Integration/Integration.jl")
include("IntegralEquations/IntegralEquations.jl")
include("Mesh/Mesh.jl")
include("ParametricEntities/ParametricEntities.jl")
include("Trees/Trees.jl")
include("IO/IO.jl")

function __init__()
    # handling of optional dependencies
    @require WriteVTK="64499a7a-5c06-52f2-abe2-ccb03c286192" include("IO/vtkIO.jl")
    @require Gmsh = "705231aa-382f-11e9-3f0c-b7cb4346fdeb" include("IO/gmshIO.jl")
end

end
