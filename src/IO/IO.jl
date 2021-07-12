module IO

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

include("plotsIO.jl")

function __init__()
    # if GmshSDK is available, load gmshIO
    try
        @eval using GmshSDK
        include("scr/IO/gmshIO.jl")
        @info "GmshSDK found in current environment. Loading `IO/gmshIO.jl`"
    catch
        @info "GmshSDK not found in current environemnt. Related funtionality will be unavailable."
    end
    # if WriteVTK is available, load vtkIO
    try
        @eval using WriteVTK
        include("src/IO/vtkIO.jl")
        @info "WriteVTK found in current environment. Loading `IO/vtkIO.jl`"
    catch
        @info "WriteVTK not found in current environemnt. Related funtionality will be unavailable."
    end
end

end # module
