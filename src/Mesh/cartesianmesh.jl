"""
    struct UniformCartesianMesh{N,T} <: AbstractMesh{N,T}

An `N`-dimensional cartesian grid given as the tensor-product of `N`
one-dimensional `LinRange{T}` grids. 

Iterating over a `UniformCartesianMesh` generates the elements which compose the mesh;
i.e. the `HyperRectangle` cells.
"""
struct UniformCartesianMesh{N,T} <: AbstractMesh{N,T}
    domain::HyperRectangle{N,T}  # stores the `low_corner` and the `high_corner` of the `UniformCartesianMesh`
    sz::NTuple{N,Int}            # number of `HyperRectangle` cells per dimension
end

Interpolation.low_corner(g::UniformCartesianMesh)  = low_corner(g.domain)
Interpolation.high_corner(g::UniformCartesianMesh) = high_corner(g.domain)

Base.size(g::UniformCartesianMesh)   = g.sz
Base.size(g::UniformCartesianMesh,i) = g.sz[i]

function grids(g::UniformCartesianMesh{N}) where N
    lc = low_corner(g)
    hc = high_corner(g)
    len = size(g) .+ 1   # need n+1 points for n cells
    return ntuple(N) do i
        LinRange(lc[i],hc[i],len[i])
    end
end
function grids(g::UniformCartesianMesh,dim)
    start = low_corner(g)[dim]
    stop = high_corner(g)[dim]
    len = size(g,dim) + 1   # need n+1 points for n cells
    return LinRange(start,stop,len)
end

Base.keys(::UniformCartesianMesh{N,T}) where {N,T} = (HyperRectangle{N,T},)

function Base.step(m::UniformCartesianMesh{N}) where N
    lc = low_corner(m)
    hc = high_corner(m)
    sz = size(m)
    return ntuple(N) do i 
        (hc[i]-lc[i])/sz[i]
    end     
end

"""
    UniformCartesianMesh(domain::HyperRectangle,sz::NTuple)
    UniformCartesianMesh(domain::HyperRectangle;step::NTuple)

Construct a uniform `UniformCartesianMesh` with `sz[d]` elements along dimension
`d`. If the kwarg `step` is passed, construct a `UniformCartesianMesh` with
elements of approximate size `step`.
"""
UniformCartesianMesh(domain::HyperRectangle{N,T},sz::Int) where {N,T} = UniformCartesianMesh(domain,ntuple(i->sz,N))

function UniformCartesianMesh(grids::NTuple{N,LinRange{T}}) where {N,T}
    lc = SVector{N,T}(minimum(r) for r in grids)
    hc = SVector{N,T}(maximum(r) for r in grids)
    sz = length.(grids) .- 1   # need n+1 points for n cells
    domain = HyperRectangle(lc,hc)
    return UniformCartesianMesh(domain,sz)
end

# in case you pass arguments like UniformCartesianMesh(xgrid,ygrid), convert
# them to a Tuple
UniformCartesianMesh(grids::Vararg{LinRange{T}}) where {T} = UniformCartesianMesh(Tuple(grids))

function UniformCartesianMesh(domain::HyperRectangle{N};step::NTuple{N}) where {N}
    lc = low_corner(domain)
    hc = high_corner(domain)
    sz = ntuple(N) do i
        (hc[i] - lc[i]) / step[i] |> ceil |> Int
    end
    return UniformCartesianMesh(domain,sz)
end

ambient_dimension(::UniformCartesianMesh{N}) where {N} = N

xgrid(g::UniformCartesianMesh) = grids(g,1)
ygrid(g::UniformCartesianMesh) = grids(g,2)
zgrid(g::UniformCartesianMesh) = grids(g,3)

# implement ElementIterator interface to UniformCartesianMesh

function Base.size(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh})
    g = mesh(iter)
    return size(g)
end

Base.length(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh}) = prod(size(iter))

function Base.CartesianIndices(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh})
    CartesianIndices(size(iter))
end

function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh}, I::CartesianIndex)
    m = mesh(iter)
    N = ambient_dimension(m)
    @assert N == length(I)
    _grids = grids(m)
    lc = ntuple(N) do dim
        i = I[dim]
        _grids[dim][i]
    end
    hc = ntuple(N) do dim
        i = I[dim] + 1
        _grids[dim][i]
    end
    return HyperRectangle(lc, hc)
end
function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh},I...)
    return iter[CartesianIndex(I)]
end
function Base.getindex(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh}, i::Int)
    I = CartesianIndices(iter)
    return iter[I[i]]
end

function Base.iterate(iter::ElementIterator{<:HyperRectangle,<:UniformCartesianMesh},state=1)
    state > length(iter) && (return nothing)
    return iter[state], state + 1
end

# since UniformCartesianMesh has only one elment type, calling ElementIterator without
# specifying the has clear sense
function ElementIterator(m::UniformCartesianMesh)
    E = keys(m) |> first
    return ElementIterator(m,E)
end

# NodeIterator over UniformCartesianMesh

function Base.IteratorSize(::Type{NodeIterator{UniformCartesianMesh{N,T}}}) where {N,T}
    return Base.HasShape{N}()
end

function Base.size(iter::NodeIterator{<:UniformCartesianMesh})
    g = mesh(iter)
    return size(g) .+ 1   # need n+1 points for n cells
end

Base.length(iter::NodeIterator{<:UniformCartesianMesh}) = prod(size(iter))

function Base.CartesianIndices(iter::NodeIterator{<:UniformCartesianMesh})
    return CartesianIndices(size(iter))
end

function Base.getindex(iter::NodeIterator{<:UniformCartesianMesh}, I::CartesianIndex)
    m = mesh(iter)
    N = ambient_dimension(m)
    @assert N == length(I)
    _grids = grids(m)
    return svector(N) do dim
        i = I[dim]
        _grids[dim][i]
    end
end
function Base.getindex(iter::NodeIterator{<:UniformCartesianMesh},I...)
    return iter[CartesianIndex(I)]
end
function Base.getindex(iter::NodeIterator{<:UniformCartesianMesh}, i::Int)
    I = CartesianIndices(iter)
    return iter[I[i]]
end

function Base.iterate(iter::NodeIterator{<:UniformCartesianMesh},state=1)
    state > length(iter) && (return nothing)
    return iter[state], state + 1
end
