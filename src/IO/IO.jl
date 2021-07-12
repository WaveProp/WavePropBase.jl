module IO

using StaticArrays
using Printf
using RecipesBase
using OrderedCollections
using Requires # for conditional loading

using WavePropBase
using WavePropBase.Utils
using WavePropBase.Geometry
using WavePropBase.Integration
using WavePropBase.Interpolation
using WavePropBase.Mesh

WavePropBase.@import_interface

include("plotsIO.jl")

function __init__()
    # # if GmshSDK is available, load gmshIO
    # @require GmshSDK include("scr/IO/gmshIO.jl")
    # if WriteVTK is available, load vtkIO
    @require WriteVTK="64499a7a-5c06-52f2-abe2-ccb03c286192" include("vtkIO.jl")
end

end # module
