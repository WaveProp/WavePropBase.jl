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

Base.size(g::CartesianMesh)   = length.(g.grids) .- 1

Base.length(g::CartesianMesh) = prod(size(g))

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
# iterate over. Below we define `getindex` for a `CartesianIndex`, which
# instantiates a `HyperRectangle`, and use that as the base implementation for
# the other `getindex` and iterate

Base.CartesianIndices(g::CartesianMesh) = CartesianIndices(size(g))

function Base.getindex(g::CartesianMesh, I::CartesianIndex)
    N = ambient_dimension(g)
    @assert N == length(I)
    lc = ntuple(N) do dim
        i = I[dim]
        g.grids[dim][i]
    end
    hc = ntuple(N) do dim
        i = I[dim] + 1
        g.grids[dim][i]
    end
    HyperRectangle(lc, hc)
end

Base.getindex(g::CartesianMesh,I...) = g[CartesianIndex(I)]
function Base.getindex(g::CartesianMesh, i::Int)
    I = CartesianIndices(g)
    g[I[i]]
end

function Base.iterate(m::CartesianMesh, state=1)
    state > length(m) && (return nothing)
    m[state], state + 1
end
