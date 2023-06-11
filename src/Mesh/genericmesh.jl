"""
    struct GenericMesh{N,T} <: AbstractMesh{N,T}

Data structure representing a generic mesh in an ambient space of dimension `N`,
with data of type `T`.

The `GenericMesh` can, in principle, store elements of any type. Those are given
as a key in the `elements` dictionary, and the value is a data structure which
is capable of reconstructing the elements. For example, for a Lagrange element
of described by `p` nodes, the value is a `p×Nel` matrix of integer, where each
columns is a list of tags for the nodes of the element. The nodes are stored in
the `nodes` field.
"""
Base.@kwdef struct GenericMesh{N,T} <: AbstractMesh{N,T}
    nodes::Vector{SVector{N,T}} = Vector{SVector{N,T}}()
    # for each element type (key), store the data required to reconstruct the
    # elements (value). The nature of the value depends on the element type.
    etype2data::Dict{DataType,Any} = Dict{DataType,Any}()
    # mapping from entity to a dict containing (etype=>tags)
    ent2tagsdict::Dict{AbstractEntity,Dict{DataType,Vector{Int}}} = Dict{AbstractEntity,
                                                                         Dict{DataType,
                                                                              Vector{Int}}}()
end

nodes(m::GenericMesh) = m.nodes

etype2data(m::GenericMesh) = m.etype2data

ent2tagsdict(m::GenericMesh) = m.ent2tagsdict

Base.keys(m::GenericMesh) = keys(etype2data(m))

entities(m::GenericMesh) = collect(keys(ent2tagsdict(m)))

domain(m::GenericMesh) = Domain(entities(m))

# implement the interface for ElementIterator of lagrange elements on a generic
# mesh. The elements are constructed on the flight based on the global nodes and
# the connectivity list stored
function Base.size(iter::ElementIterator{<:LagrangeElement,<:GenericMesh})
    msh = mesh(iter)
    E   = eltype(iter)
    tags = msh.etype2data[E]::Matrix{Int}
    _, Nel = size(tags)
    return (Nel,)
end

function Base.getindex(iter::ElementIterator{<:LagrangeElement,<:GenericMesh}, i::Int)
    E = eltype(iter)
    M = mesh(iter)
    tags = M.etype2data[E]::Matrix{Int}
    node_tags = view(tags, :, i)
    vtx = view(M.nodes, node_tags)
    el = E(vtx)
    return el
end

function Base.iterate(iter::ElementIterator{<:LagrangeElement,<:GenericMesh}, state=1)
    state > length(iter) && (return nothing)
    return iter[state], state + 1
end

# element iterator interface for parametric elements. The elements are stored
# directly in the etype2data field
function Base.size(iter::ElementIterator{<:ParametricElement,<:GenericMesh})
    E = eltype(iter)
    M = mesh(iter)
    els = M.etype2data[E]::Vector{E}
    return (length(els),)
end

function Base.getindex(iter::ElementIterator{<:ParametricElement,<:GenericMesh}, i::Int)
    E = eltype(iter)
    M = mesh(iter)
    els = M.etype2data[E]::Vector{E}
    return els[i]
end

function Base.iterate(iter::ElementIterator{<:ParametricElement,<:GenericMesh}, state=1)
    state > length(iter) && (return nothing)
    return iter[state], state + 1
end

"""
    dom2elt(m::GenericMesh,Ω,E)

Compute the element indices `idxs` of the elements of type `E` composing `Ω`, so
that `m[E][idxs]` gives all the elements of type `E` meshing `Ω`.
"""
function dom2elt(m::GenericMesh, Ω, E::DataType)
    idxs = Int[]
    for ent in entities(Ω)
        tags = get(m.ent2tagsdict[ent], E, Int[])
        append!(idxs, tags)
    end
    return idxs
end

"""
    dom2elt_dict(m::GenericMesh,Ω)

Return a `Dict` with keys being the element types of `m`, and values being the
indices of the elements in `Ω` of that type.

See also: [`dom2elt`](@ref)
"""
function dom2elt_dict(m::GenericMesh, Ω)
    dict = Dict{DataType,Vector{Int}}()
    for E in keys(m)
        tags = dom2elt(m, Ω, E)
        if !isempty(tags)
            dict[E] = tags
        end
    end
    return dict
end

# convert a mesh to 2d by ignoring third component. Note that this also requires
# converting various element types to their 2d counterpart. These are needed
# because some meshers like gmsh always create three-dimensional objects, so we
# must convert after importing the mesh
function convert_to_2d(mesh::GenericMesh{3})
    @assert all(E -> geometric_dimension(domain(E)) < 3, keys(mesh))
    T = primitive_type(mesh)
    # create new dictionaries for elements and ent2tagsdict with 2d elements as keys
    els = empty(etype2data(mesh))
    e2t = empty(ent2tagsdict(mesh))
    for (E, tags) in etype2data(mesh)
        E2d = convert_to_2d(E)
        els[E2d] = tags
    end
    for (ent, dict) in ent2tagsdict(mesh)
        new_dict = empty(dict)
        for (E, tags) in dict
            E2d = convert_to_2d(E)
            new_dict[E2d] = tags
        end
        e2t[ent] = new_dict
    end
    # construct new 2d mesh
    return GenericMesh{2,T}(;
                            nodes=[x[1:2] for x in nodes(mesh)],
                            etype2data=els,
                            ent2tagsdict=e2t)
end

function convert_to_2d(::Type{LagrangeElement{R,N,SVector{3,T}}}) where {R,N,T}
    return LagrangeElement{R,N,SVector{2,T}}
end
