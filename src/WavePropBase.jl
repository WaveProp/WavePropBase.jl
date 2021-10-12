module WavePropBase

include("interface.jl")

include("Utils/Utils.jl")

include("Geometry/Geometry.jl")

include("Interpolation/Interpolation.jl")

include("Integration/Integration.jl")

include("Mesh/Mesh.jl")

include("Trees/Trees.jl")

include("IO/IO.jl")

@export_interface

# export modules. IO module is not exported due to name conflict.
export
    Utils,
    Geometry,
    Trees,
    Interpolation,
    Integration,
    Mesh

end
