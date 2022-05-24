"""
    abstract type AbstractReferenceShape{N}

A reference domain/shape in `ℜᴺ`.

Used mostly for defining more complex shapes as transformations mapping an
`AbstractReferenceShape` into some region of `ℜᴹ`.

See e.g. [`ReferenceLine`](@ref) or [`ReferenceTriangle`](@ref) for some
examples of concrete subtypes.
"""
abstract type AbstractReferenceShape{N} end

"""
    ambient_dimension(x)

Dimension of the ambient space where `x` lives. For geometrical objects this can
differ from its [`geometric_dimension`](@ref); for example a triangle in `ℝ³` has
ambient dimension `3` but geometric dimension `2`, while a curve in `ℝ³` has
ambient dimension 3 but geometric dimension 1.
"""
ambient_dimension(::SType{<:AbstractReferenceShape{N}}) where {N}    = N


geometric_dimension(::SType{<:AbstractReferenceShape{N}}) where {N}  = N


dimension(::SType{<:AbstractReferenceShape{N}}) where {N}  = N

"""
    struct ReferencePoint

Singleton type representing a reference zero-dimensional entitty (i.e. a point).
"""
struct ReferencePoint <: AbstractReferenceShape{0}
end
number_of_nodes(::SType{ReferencePoint}) = 1

"""
    struct ReferenceLine

Singleton type representing the `[0,1]` segment.
"""
struct ReferenceLine <: AbstractReferenceShape{1}
end
Base.in(x,::ReferenceLine)    = 0 ≤ x[1] ≤ 1
center(::Type{ReferenceLine}) = 0.5
center(::ReferenceLine)       = 0.5
number_of_nodes(::Type{ReferenceLine}) = 2
number_of_nodes(::ReferenceLine)       = 2

vertices(ln::ReferenceLine) = SVector(0), SVector(1)

"""
    struct ReferenceTriangle

Singleton type representing the triangle with vertices `(0,0),(0,1),(1,0)`
"""
struct ReferenceTriangle <: AbstractReferenceShape{2}
end
Base.in(x,::ReferenceTriangle) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1 - x[1]
number_of_nodes(::Type{ReferenceTriangle}) = 3
number_of_nodes(::ReferenceTriangle)       = 3

vertices(::ReferenceTriangle) = SVector(0,0), SVector(1,0), SVector(0,1)

"""
    struct ReferenceHyperCube{N}

Singleton type representing the axis-aligned hypercube in `N` dimensions with
the lower corner at the origin and the upper-corner at `(1,1,…,1)`.
"""
struct ReferenceHyperCube{N} <: AbstractReferenceShape{N}
end
Base.in(x,::ReferenceHyperCube)    = all(0 ≤ xi ≤ 1 for xi in x)

center(::SType{ReferenceHyperCube{N}}) where {N} = svector(i->0.5,N)

number_of_nodes(::SType{ReferenceHyperCube{N}}) where {N} = 2^N

"""
    const ReferenceSquare = ReferenceHyperCube{2}

Singleton type representing the square with vertices `(0,0),(0,1),(1,1),(1,0)`
"""
const ReferenceSquare = ReferenceHyperCube{2}

vertices(sq::ReferenceSquare) = SVector(0,0), SVector(1,0), SVector(1,1), SVector(0,1)

"""
    struct ReferenceTetrahedron

Singleton type representing the tetrahedron with vertices `(0,0,0),(0,0,1),(0,1,0),(1,0,0)`
"""
struct ReferenceTetrahedron <: AbstractReferenceShape{3}
end
Base.in(x,::ReferenceTetrahedron) = 0 ≤ x[1] ≤ 1 &&
                                    0 ≤ x[2] ≤ 1 - x[1] &&
                                    0 ≤ x[3] ≤ 1 - x[1] - x[2]
number_of_nodes(::Type{ReferenceTetrahedron}) = 4
number_of_nodes(::ReferenceTetrahedron)       = 4

# TODO: generalize structs above to `ReferenceSimplex{N}` and
# `ReferenceCuboid{N}`
