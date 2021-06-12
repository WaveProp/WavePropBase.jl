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
CartesianMesh(args...) = CartesianMesh(args)

grids(g::CartesianMesh)     = g.grids
grids(g::CartesianMesh,dim) = g.grids[dim]

Base.keys(m::CartesianMesh{N,T}) where {N,T} = (HyperRectangle{N,T},)

"""
    CartesianMesh(;domain::HyperRectangle,sz::NTuple)

Construct a uniform `CartesianMesh` with `sz[d]` elements along dimension `d`.
"""
function CartesianMesh(domain::HyperRectangle{N,T},sz::NTuple{N}) where {N,T}
    lc = low_corner(domain)
    uc = high_corner(domain)
    grids1d = ntuple(N) do n
        # to have sz elements, need sz+1 points in LinRange
        npts = sz[n] + 1
        LinRange{T}(lc[n], uc[n], npts) |> collect
    end
    CartesianMesh(grids1d)
end
CartesianMesh(domain::HyperRectangle{N,T},sz::Int) where {N,T} = CartesianMesh(domain,ntuple(i->sz,N))
CartesianMesh(;domain,sz) = CartesianMesh(domain,sz)

ambient_dimension(g::CartesianMesh{N}) where {N} = N

xgrid(g::CartesianMesh) = g.grids[1]

ygrid(g::CartesianMesh) = g.grids[2]

zgrid(g::CartesianMesh) = g.grids[3]

# implement ElementIterator interface to CartesianMesh

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

# since CartesianMesh has only one elment type, calling ElementIterator without
# specifying the has clear sense
function ElementIterator(m::CartesianMesh)
    E = keys(m) |> first
    ElementIterator(m,E)
end

# NodeIterator over CartesianMesh

function Base.IteratorSize(iter::Type{NodeIterator{CartesianMesh{N,T}}}) where {N,T}
    Base.HasShape{N}()
end

function Base.size(iter::NodeIterator{<:CartesianMesh})
    g = mesh(iter)
    length.(g.grids)
end

Base.length(iter::NodeIterator{<:CartesianMesh}) = prod(size(iter))

function Base.CartesianIndices(iter::NodeIterator{<:CartesianMesh})
    CartesianIndices(size(iter))
end

function Base.getindex(iter::NodeIterator{<:CartesianMesh}, I::CartesianIndex)
    m = mesh(iter)
    N = ambient_dimension(m)
    @assert N == length(I)
    svector(N) do dim
        i = I[dim]
        m.grids[dim][i]
    end
end
function Base.getindex(iter::NodeIterator{<:CartesianMesh},I...)
    iter[CartesianIndex(I)]
end
function Base.getindex(iter::NodeIterator{<:CartesianMesh}, i::Int)
    I = CartesianIndices(iter)
    iter[I[i]]
end

function Base.iterate(iter::NodeIterator{<:CartesianMesh},state=1)
    state > length(iter) && (return nothing)
    iter[state], state + 1
end
