module WavePropBase

using StaticArrays

include("referencedomain.jl")

export
    # types
    ReferenceLine,
    ReferenceTriangle,
    ReferenceSquare,
    ReferenceTetrahedron,
    # methods
    ambient_dimension,
    geometric_dimension

end
