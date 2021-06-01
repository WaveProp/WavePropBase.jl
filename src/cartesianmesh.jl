"""
    struct CartesianMesh{N,T} <: AbstractMesh{N,T}

An `N`-dimensional cartesian grid given as the tensor-product of `N`
one-dimensional `Vector{T}` grids.

Iterating over a `CartesianMesh` generates the elements which compose the mesh;
i.e. the `HyperRectangle` cells.
"""
struct CartesianMesh{N,T} <: AbstractMesh{N,T}
    grids::NTuple{N,Vector{T}}
end

grids(g::CartesianMesh)     = g.grids
grids(g::CartesianMesh,dim) = g.grids[dim]

Base.keys(m::CartesianMesh{N,T}) where {N,T} = (HyperRectangle{N,T},)

"""
    CartesianMesh(;domain::HyperRectangle,sz::NTuple)

Construct a uniform `CartesianMesh` with `sz[d]` elements along dimension `d`.
"""
function CartesianMesh(;domain::HyperRectangle{N,T},sz::NTuple{N}) where {N,T}
    lc = low_corner(domain)
    uc = high_corner(domain)
    grids1d = ntuple(N) do n
        # to have sz elements, need sz+1 points in LinRange
        npts = sz[n] + 1
        LinRange{T}(lc[n], uc[n], npts) |> collect
    end
    CartesianMesh(grids1d)
end

ambient_dimension(g::CartesianMesh{N}) where {N} = N

xgrid(g::CartesianMesh) = g.grids[1]

ygrid(g::CartesianMesh) = g.grids[2]

zgrid(g::CartesianMesh) = g.grids[3]

# the cartesian structure of `CartesianMesh` makes it particularly easy to
# iterate over using CartesianIndices. Iterating over linear indices then simply
# convert to a CartesianIndex.

function Base.size(iter::ElementIterator{<:HyperRectangle,<:CartesianMesh})
    g = mesh(iter)
    length.(g.grids) .- 1
end

Base.length(iter::ElementIterator{<:HyperRectangle,<:CartesianMesh}) = prod(size(iter))

function Base.CartesianIndices(iter::ElementIterator{<:HyperRectangle,<:CartesianMesh})
    CartesianIndices(size(iter))
end

function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:CartesianMesh}, I::CartesianIndex)
    m = mesh(iter)
    N = ambient_dimension(m)
    @assert N == length(I)
    lc = ntuple(N) do dim
        i = I[dim]
        m.grids[dim][i]
    end
    hc = ntuple(N) do dim
        i = I[dim] + 1
        m.grids[dim][i]
    end
    HyperRectangle(lc, hc)
end
function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:CartesianMesh},I...)
    iter[CartesianIndex(I)]
end
function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:CartesianMesh}, i::Int)
    I = CartesianIndices(iter)
    iter[I[i]]
end

function Base.iterate(iter::ElementIterator{<:HyperRectangle,<:CartesianMesh},state=1)
    state > length(iter) && (return nothing)
    iter[state], state + 1
end
