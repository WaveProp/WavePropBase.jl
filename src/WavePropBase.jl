module WavePropBase

include("interface.jl")

include("Utils/Utils.jl")

include("Geometry/Geometry.jl")

include("Trees/Trees.jl")

include("Interpolation/Interpolation.jl")

include("Integration/Integration.jl")

include("Mesh/Mesh.jl")

include("IO/IO.jl")

include("Simulation/Simulation.jl")

@export_interface

# export modules. IO module is not exported due to name conflict.
export
    Utils,
    Geometry,
    Interpolation,
    Integration,
    Mesh,
    Simulation

end
