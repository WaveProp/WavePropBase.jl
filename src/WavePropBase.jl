module WavePropBase

const PROJECT_ROOT =  pkgdir(WavePropBase)

using StaticArrays
using LinearAlgebra
using AbstractTrees
using RecipesBase
using Requires # for conditional loading of vtkIO
using Printf
using Statistics: median

include("Utils/Utils.jl")
include("Geometry/Geometry.jl")
include("Interpolation/Interpolation.jl")
include("Integration/Integration.jl")
include("Mesh/Mesh.jl")
include("Trees/Trees.jl")
include("IO/IO.jl")

function __init__()
    # if WriteVTK is available, include vtkIO
    @require WriteVTK="64499a7a-5c06-52f2-abe2-ccb03c286192" begin
        @info "including vtkIO.jl from WavePropBase/IO"
        include("IO/vtkIO.jl")
    end
end

end
