module WavePropBase

using LinearAlgebra
using StaticArrays

include("utils.jl")
include("interface.jl")
include("referenceshapes.jl")
include("hyperrectangle.jl")
include("entities.jl")
include("domain.jl")
include("lagrangeinterp.jl")
include("element.jl")
include("abstractmesh.jl")
include("cartesianmesh.jl")
include("genericmesh.jl")

export
    # Abstract types
    AbstractMesh,
    AbstractEntity,
    ReferenceShape,
    # Concrete types
    HyperRectangle,
    ReferenceLine,
    ReferenceTriangle,
    ReferenceSquare,
    Domain,
    CartesianMesh,
    # Type aliases
    Point1D,
    Point2D,
    Point3D

end
