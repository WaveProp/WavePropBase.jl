"""
    abstract type ReferenceShape{N}

Singleton types defining a referene domain/shape in `ℜᴺ`. Typically used to define
more complex shapes as transformations mapping an `ReferenceShape` into
some region of `ℜᴹ` through a given map.

See e.g. [`ReferenceLine`](@ref) or [`ReferenceTriangle`](@ref) for
examples of concrete subtypes.
"""
abstract type ReferenceShape{N} end

ambient_dimension(::SType{<:ReferenceShape{N}}) where {N}    = N
geometric_dimension(::SType{<:ReferenceShape{N}}) where {N}  = N

"""
    struct ReferenceLine

Singleton type representing the `[0,1]` segment.
"""
struct ReferenceLine <: ReferenceShape{1}
end
Base.in(x,::ReferenceLine)    = 0 ≤ x[1] ≤ 1
getcenter(::SType{ReferenceLine}) = 0.5

getvertices(::SType{ReferenceLine}) = SVector(0), SVector(1)

"""
    struct ReferenceTriangle

Singleton type representing the triangle with vertices `(0,0),(0,1),(1,0)`
"""
struct ReferenceTriangle <: ReferenceShape{2}
end
Base.in(x,::ReferenceTriangle) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1 - x[1]

getvertices(::SType{ReferenceTriangle}) = SVector(0,0), SVector(1,0), SVector(0,1)

"""
    struct ReferenceSquare

Singleton type representing the square with vertices `(0,0),(0,1),(1,1),(1,0)`
"""
struct ReferenceSquare <: ReferenceShape{2}
end
Base.in(x,::ReferenceSquare)        = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1
getcenter(::SType{ReferenceSquare}) = SVector(0.5,0.5)

getvertices(::SType{ReferenceSquare}) = SVector(0,0), SVector(1,0), SVector(1,1), SVector(0,1)

"""
    struct ReferenceTetrahedron

Singleton type representing the tetrahedron with vertices `(0,0,0),(0,0,1),(0,1,0),(1,0,0)`
"""
struct ReferenceTetrahedron <: ReferenceShape{3}
end
Base.in(x,::ReferenceTetrahedron) = 0 ≤ x[1] ≤ 1 &&
                                    0 ≤ x[2] ≤ 1 - x[1] &&
                                    0 ≤ x[3] ≤ 1 - x[1] - x[2]
