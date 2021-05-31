"""
    abstract type AbstractMesh{N,T}

An abstract mesh structure in dimension `N` with primite data of type `T` (e.g.
`Float64` for double precision representation).

The `AbstractMesh` interface expects the mesh object to overload `Base.keys(m)`
to return all the element types composing the mesh. Use a key to index the mesh
is expected to return an iterator over all elements of that type. Finally,
calling `mesh[E][i]` should return the *i-the* element of type `E`.
"""
abstract type AbstractMesh{N,T} end

ambient_dimension(M::AbstractMesh{N}) where {N} = N

# we define the geometric dimension of the mesh to be the largest of the geometric
# dimension of its entities.
geometric_dimension(M::AbstractMesh) = maximum(x -> geometric_dimension(x), etypes(M))

Base.eltype(M::AbstractMesh{N,T}) where {N,T} = T

Base.length(mesh::AbstractMesh) = length(nodes(mesh))

"""
    struct ElementIterator{E,M}

Return an iterator for iterating over all elements of type `E` on mesh objects of
type `M`.
"""
struct ElementIterator{E,M}
    mesh::M
end

mesh(iter::ElementIterator) = iter.mesh

Base.eltype(::SType{ElementIterator{E}}) where {E} = E

ElementIterator{E}(mesh::M) where {E,M <: AbstractMesh} = ElementIterator{E,M}(mesh)

ElementIterator(mesh,E) = ElementIterator{E}(mesh)

Base.iterate(iter::ElementIterator, state=1) = iterate(mesh(iter),state)

Base.size(iter::ElementIterator{<:Any,<:CartesianMesh})   = size(mesh(iter))

Base.length(iter::ElementIterator{<:Any,<:CartesianMesh}) = length(mesh(iter))
