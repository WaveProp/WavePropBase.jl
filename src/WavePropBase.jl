module WavePropBase

using LinearAlgebra
using StaticArrays
using RecipesBase

include("interface.jl")

include("Utils/Utils.jl")

include("Geometry/Geometry.jl")

include("Interpolation/Interpolation.jl")

include("Integration/Integration.jl")

include("Mesh/Mesh.jl")

include("IO/IO.jl")

@export_interface

end
