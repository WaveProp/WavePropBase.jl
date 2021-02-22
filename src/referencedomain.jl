"""
    abstract type AbstractReferenceDomain{N}

A reference domain/shape in `ℜᴺ`.

Used mostly for defining more complex shapes as transformations mapping an
`AbstractReferenceDomain` into some region of `ℜᴹ`, where `M` need not be the
same as `N`.

See e.g. [`ReferenceLine`](@ref) or [`ReferenceTriangle`](@ref) for some
examples of concrete subtypes.
"""
abstract type AbstractReferenceDomain{N} end

ambient_dimension(::Type{<:AbstractReferenceDomain{N}}) where {N}   = N
geometric_dimension(::Type{<:AbstractReferenceDomain{N}}) where {N} = N
ambient_dimension(::AbstractReferenceDomain{N}) where {N}           = N
geometric_dimension(::AbstractReferenceDomain{N}) where {N}         = N

"""
    struct ReferenceLine <: AbstractReferenceDomain{1}

Singleton type representing the `[0,1]` segment.
"""
struct ReferenceLine <: AbstractReferenceDomain{1}
end
Base.in(x,::ReferenceLine)    = 0 ≤ x[1] ≤ 1

"""
    struct ReferenceTriangle <: AbstractReferenceDomain{2}

Singleton type representing the triangle with vertices `(0,0),(0,1),(1,0)`
"""
struct ReferenceTriangle <: AbstractReferenceDomain{2}
end
Base.in(x,::ReferenceTriangle) = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1 - x[1]

"""
    struct ReferenceSquare <: AbstractReferenceDomain{2}

Singleton type representing the square with vertices `(0,0),(0,1),(1,1),(1,0)`
"""
struct ReferenceSquare <: AbstractReferenceDomain{2}
end
Base.in(x,::ReferenceSquare)    = 0 ≤ x[1] ≤ 1 && 0 ≤ x[2] ≤ 1

"""
    struct ReferenceTetrahedron <: AbstractReferenceDomain{3}

Singleton type representing the tetrahedron with vertices `(0,0,0),(0,0,1),(0,1,0),(1,0,0)`
"""
struct ReferenceTetrahedron <: AbstractReferenceDomain{3}
end
Base.in(x,::ReferenceTetrahedron) = 0 ≤ x[1] ≤ 1 &&
                                    0 ≤ x[2] ≤ 1 - x[1] &&
                                    0 ≤ x[3] ≤ 1 - x[1] - x[2]
