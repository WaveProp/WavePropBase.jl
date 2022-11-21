"""
    abstract type AbstractMesh{N,T}

An abstract mesh structure in dimension `N` with primite data of type `T` (e.g.
`Float64` for double precision representation).

Concrete subtypes of `AbstractMesh` should implement [`ElementIterator`](@ref)
for accessing the mesh elements, and/or [`NodeIterator`](@ref) for accesing the
mesh nodes.

# See also: [`GenericMesh`](@ref), [`UniformCartesianMesh`](@ref)
"""
abstract type AbstractMesh{N,T} end

ambient_dimension(M::AbstractMesh{N}) where {N} = N

# we define the geometric dimension of the parent to be the largest of the geometric
# dimension of its entities.
geometric_dimension(M::AbstractMesh) = maximum(x -> geometric_dimension(x), entities(M))

primitive_type(M::AbstractMesh{N,T}) where {N,T} = T

function Base.show(io::IO, msh::AbstractMesh)
    print(io, "$(typeof(msh)) containing:")
    for E in keys(msh)
        iter = ElementIterator(msh, E)
        print(io, "\n\t $(length(iter)) elements of type ", E)
    end
    return io
end

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

ElementIterator{E}(parent::M) where {E,M<:AbstractMesh} = ElementIterator{E,M}(parent)

ElementIterator(parent, E) = ElementIterator{E}(parent)

# indexing a mesh with an element type return an interator over that element
# (like a dict)
Base.getindex(m::AbstractMesh, E::DataType) = ElementIterator(m, E)

# TODO: rename this to elements and warn against the fact that it returns a
# possibly type unstable vector
function all_elements(msh::AbstractMesh)
    T = Union{keys(msh)...}
    els = Vector{T}()
    for E in keys(msh)
        for el in msh[E]
            push!(els, el)
        end
    end
    return els
end

function Base.length(iter::ElementIterator)
    return prod(size(iter))
end

"""
    struct NodeIterator{M}

Iterator for all the nodes in a mesh of type `M`.
"""
struct NodeIterator{M}
    mesh::M
end

mesh(iter::NodeIterator) = iter.mesh

# nodes of `AbstractMesh{N,T}` should be points of type `SVector{N,T}`
Base.eltype(::SType{NodeIterator{<:AbstractMesh{N,T}}}) where {N,T} = SVector{N,T}

"""
    near_interaction_list(X,Y::AbstractMesh; tol)

For each element `el` of type `E` in `Y`, return the indices of the points in `X` which
are closer than `tol` to the `center` of `el`.

This function returns a dictionary where e.g. `Dict[E][5] --> Vector{Int}` gives
the indices of points in `X` which are closer than `tol` to the center of the
fifth element of type `E`.
"""
function near_interaction_list(X, Y::AbstractMesh{N}; tol) where {N}
    # bounding box for target elements with a layer of width `tol` so that
    # points which are not in the bounding box but are farther than `tol` from
    # any target point
    @assert tol > 0
    bbox = HyperRectangle(X)
    bbox = HyperRectangle(low_corner(bbox).-tol, high_corner(bbox).+tol)
    # a cartesian mesh to sort the points
    msh  = UniformCartesianMesh(bbox;step=tol)
    sz   = size(msh)
    targets_in_element = sort_in_cartesian_mesh(X,msh)
    # for each element type, build the list of targets close to a given element
    dict = Dict{DataType,Vector{Vector{Int}}}()
    for E in keys(Y)
        nel = length(Y[E])
        dict[E] = idxs = [Int[] for _ in 1:nel]
        # TODO: function barrier needed
        for (n,el) in enumerate(Y[E])
            y = center(el)
            I = element_index_for_point(y,msh)
            isnothing(I) && continue
            It = Tuple(I)
            # index of neighboring boxes
            neighbors = ntuple(N) do d
                max(It[d]-1,1):min(It[d]+1,sz[d])
            end |> CartesianIndices
            # append targets to idxs of current element
            for J in neighbors
                haskey(targets_in_element,J) || continue
                append!(idxs[n],targets_in_element[J])
            end
        end
    end
    return dict
end

# function near_interaction_list(X, M::AbstractMesh; tol)
#     dict = Dict{DataType,Vector{Vector{Int}}}()
#     for E in keys(M)
#         pts = [center(el) for el in M[E]]
#         kdtree = KDTree(pts)
#         push!(dict,E=>inrange(kdtree, X, tol))
#     end
#     return dict
# end
