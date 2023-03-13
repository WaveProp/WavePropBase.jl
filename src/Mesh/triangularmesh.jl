"""
    struct TriangularMesh{N,T} <: AbstractMesh{N,T}

A simple mesh structure containing only (flat) triangles in `N` dimensions. The
underlying data type is given by `T`, which defaults to `Float64`.

`TriangularMesh`es are presented as a collection of `nodes` and a `connectivity`
matrix of size `3 × Nt`, where `Nt` is the number of triangles in the mesh.
"""
struct TriangularMesh{N,T} <: AbstractMesh{N,T}
    nodes::Vector{SVector{N,T}}
    connectivity::Matrix{Int}
end

function TriangularMesh(msh::AbstractMesh{N,T}) where {N,T}
    nodes = SVector{N,T}[]
    connectivity = Int[]
    old2new = Dict{Int,Int}() # map from old node index to new node index
    for E in keys(msh)
        if E != LagrangeTriangle{3,SVector{N,T}}
            @warn "skypping elements of type $E"
            continue
        end
        el2idxs = msh.elements[E] # 3 × n connectivity matrix
        for n in 1:size(el2idxs, 2)
            for i in 1:3
                idx = el2idxs[i, n]
                if idx ∉ keys(old2new) # point is new
                    push!(nodes, msh.nodes[idx])
                    old2new[idx] = length(nodes)
                end
                push!(connectivity, old2new[idx])
            end
        end
    end
    connectivity = reshape(connectivity, 3, :)
    return TriangularMesh{N,T}(nodes, connectivity)
end

Base.keys(::TriangularMesh{N,T}) where {N,T} = (LagrangeTriangle{3,SVector{N,T}},)

function Base.size(iter::ElementIterator{<:LagrangeTriangle{3,SVector{N,T}},
                                         TriangularMesh{N,T}}) where {N,T}
    msh = mesh(iter)
    tags = msh.connectivity
    _, Nel = size(tags)
    return (Nel,)
end

function Base.getindex(iter::ElementIterator{<:LagrangeTriangle{3,SVector{N,T}},
                                             TriangularMesh{N,T}}, i::Int) where {N,T}
    msh = mesh(iter)
    tags = msh.connectivity
    node_tags = view(tags, :, i)
    vtx = view(msh.nodes, node_tags)
    E = eltype(iter)
    el = E(vtx)
    return el
end

function Base.iterate(iter::ElementIterator{<:LagrangeTriangle,<:TriangularMesh}, state=1)
    state > length(iter) && (return nothing)
    return iter[state], state + 1
end
