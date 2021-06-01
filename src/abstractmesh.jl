"""
    abstract type AbstractMesh{N,T}

An abstract parent structure in dimension `N` with primite data of type `T` (e.g.
`Float64` for double precision representation).

The `AbstractMesh` interface expects the following methods to be implemented:

- `keys(msh)` : return a list of the element types composing the mesh.
- `msh[E]`    : return an [`ElementIterator`](@ref) for the mesh elements of type `E`

# See also: [`ElementIterator`](@ref)
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

Return an iterator for all elements of type `E` on a mesh of type `M`.

Besides the methods listed in the [iterator
iterface](https://docs.julialang.org/en/v1/manual/interfaces/) of `Julia`, some
functions also require the `getindex(iter,i::Int)` method for accessing the
`i`-th element directly.
"""
struct ElementIterator{E,M}
    mesh::M
end

mesh(iter::ElementIterator) = iter.mesh

Base.eltype(::SType{ElementIterator{E}}) where {E} = E

ElementIterator{E}(parent::M) where {E,M <: AbstractMesh} = ElementIterator{E,M}(parent)

ElementIterator(parent,E) = ElementIterator{E}(parent)

# indexing a mesh with an element type return an interator over that element
# (like a dict)
Base.getindex(m::AbstractMesh,E::DataType) = ElementIterator(m, E)
