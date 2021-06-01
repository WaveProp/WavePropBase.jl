module WavePropBase

using LinearAlgebra
using StaticArrays
using RecipesBase

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
include("submesh.jl")
include("plotrecipes.jl")

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
    ReferenceTetrahedron,
    Domain,
    CartesianMesh,
    GenericMesh,
    ElementIterator,
    ElementaryEntity,
    LagrangeElement,
    LagrangeLine,
    LagrangeTriangle,
    LagrangeTetrahedron,
    LagrangeInterp,
    # Type aliases
    Point1D,
    Point2D,
    Point3D,
    # Global consts
    ENTITIES,
    # methods
    ambient_dimension,
    geometric_dimension,
    domain,
    entities,
    external_boundary,
    internal_boundary,
    boundary,
    tag,
    skeleton,
    clear_entities!,
    radius,
    diameter,
    center,
    bounding_box
end
