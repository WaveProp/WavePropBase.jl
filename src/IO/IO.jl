module IO

using StaticArrays
using Printf
using RecipesBase
using OrderedCollections
using Requires # for conditional loading of vtkIO

using WavePropBase
using WavePropBase.Utils
using WavePropBase.Geometry
using WavePropBase.Trees
using WavePropBase.Integration
using WavePropBase.Interpolation
using WavePropBase.Mesh

WavePropBase.@import_interface

include("plotsIO.jl")

function __init__()
    # if WriteVTK is available, include vtkIO
    @require WriteVTK="64499a7a-5c06-52f2-abe2-ccb03c286192" begin
        @info "including vtkIO.jl from WavePropBase/IO"
        include("vtkIO.jl")
    end
end

end # module
