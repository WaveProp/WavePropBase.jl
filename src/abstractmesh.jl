"""
    abstract type AbstractMesh{N,T}

An abstract parent structure in dimension `N` with primite data of type `T` (e.g.
`Float64` for double precision representation).

The `AbstractMesh` interface expects the parent object to overload `Base.keys(m)`
to return all the element types composing the parent. Use a key to index the parent
is expected to return an iterator over all elements of that type. Finally,
calling `parent[E][i]` should return the *i-the* element of type `E`.
"""
abstract type AbstractMesh{N,T} end

ambient_dimension(M::AbstractMesh{N}) where {N} = N

# we define the geometric dimension of the parent to be the largest of the geometric
# dimension of its entities.
geometric_dimension(M::AbstractMesh) = maximum(x -> geometric_dimension(x), entities(M))

Base.eltype(M::AbstractMesh{N,T}) where {N,T} = T

Base.length(parent::AbstractMesh) = length(nodes(parent))

"""
    struct ElementIterator{E,M}

Return an iterator for iterating over all elements of type `E` on a `parent`
mesh of type `M`.
"""
struct ElementIterator{E,M}
    parent::M
end

parent(iter::ElementIterator) = iter.parent

Base.eltype(::SType{ElementIterator{E}}) where {E} = E

ElementIterator{E}(parent::M) where {E,M <: AbstractMesh} = ElementIterator{E,M}(parent)

ElementIterator(parent,E) = ElementIterator{E}(parent)

Base.getindex(m::AbstractMesh,E::AbstractElement) = ElementIterator(m,E)
