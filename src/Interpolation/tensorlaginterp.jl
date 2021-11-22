"""
    struct TensorLagInterp{N,Td,T}

Generic Lagrange interpolation over an `N`-dimensional tensor grid. The
implementation uses a multidimensional generalization of the barycentric
formula.

The main constructor takes an `SVector{N,Vector{Td}}` containig the `N`
one-dimensional `nodes` and an `Array{N,T}` of the function `vals` at the tensor
product grid formed by the one-dimensional `nodes`.

# Examples:
```julia
nx = 10
ny = 12
x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
f   = (x) -> cos(x[1]*x[2])
vals = [f((x,y)) for x in x, y in y]
p   = TensorLagInterp(SVector(x,y),vals)
p((0.1,0.2)) ≈ f((0.1,0.2))
```
"""
mutable struct TensorLagInterp{N,Td,T}
    vals::Array{T,N}
    nodes1d::NTuple{N,Vector{Td}}
    weights1d::NTuple{N,Vector{Td}}
end

interpolation_nodes(p::TensorLagInterp)                 = Iterators.product(p.nodes1d...)
interpolation_nodes(p::TensorLagInterp,d::Int,i::Int)   = p.nodes1d[d][i]
function interpolation_nodes(p::TensorLagInterp{N},I::CartesianIndex{N}) where {N}
    svector(N) do d
        i = I[d]
        p.nodes1d[d][i]
    end
end

weights1d(p::TensorLagInterp) = p.weights1d
vals(p::TensorLagInterp)      = p.vals

return_type(::TensorLagInterp{_,__,T}) where {_,__,T} = T

ambient_dimension(::TensorLagInterp{N}) where {N} = N

Base.eachindex(p::TensorLagInterp)  = eachindex(vals(p))
function Base.CartesianIndices(p::TensorLagInterp)
    sz = length.(p.nodes1d)
    CartesianIndices(sz)
end

# if passed the values and nodes, construct weights automatically
function TensorLagInterp(vals::Array{T,N},nodes::NTuple{N,Vector{Td}}) where {N,Td,T}
    weights = map(barycentric_lagrange_weights,nodes)
    TensorLagInterp{N,Td,T}(vals,nodes,weights)
end
TensorLagInterp(vals::Vector,nodes::Vector) = TensorLagInterp(vals,(nodes,))
TensorLagInterp(vals::Array,nodes::NTuple)  = TensorLagInterp(vals,SVector(nodes))

# barycentric lagrange interpolation
function (p::TensorLagInterp{N,Td,T})(x::SVector) where {N,Td,T}
    return bclag_interp(x, p.vals, p.nodes1d, p.weights1d, Val{N}(), 1, length(p.vals))
end
(p::TensorLagInterp{N,Td})(x) where {N,Td} =  p(SVector{N,Td}(x))
(p::TensorLagInterp{1})(x::Number) = p((x,))

# low-level function for multidimensional barycentric lagrange interpolation
# (adapted from FastChebInterp package)
@inline function bclag_interp(x::SVector{N,Td}, vals::Array{<:Any,N},
                            nodes::NTuple{N}, weights::NTuple{N},
                            ::Val{dim}, i1, len) where {N,Td,dim}
    T   = eltype(vals)
    n   = size(vals,dim)
    @inbounds xd = x[dim]
    @inbounds W  = weights[dim]
    @inbounds X  = nodes[dim]
    num = zero(T)
    den = zero(Td)
    if dim == 1
        for i in 1:n
            ci     = vals[i1+(i-1)]
            wi     = W[i]
            x_m_xi = xd - X[i]
            if x_m_xi == 0
                return ci
            end
            num += (wi/x_m_xi)*ci
            den += (wi/x_m_xi)
        end
        return num/den
    else
        Δi = len ÷ n # column-major stride of current dimension
        # recurse down on dimension
        dim′ = Val{dim-1}()
        for i in 1:n
            ci     = bclag_interp(x, vals, nodes, weights, dim′, i1+(i-1)*Δi, Δi)
            wi     = W[i]
            x_m_xi = xd - X[i]
            if x_m_xi == 0
                return ci
            end
            num += wi/x_m_xi*ci
            den += wi/x_m_xi
        end
        return num/den
    end
end

# algorithm on https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function barycentric_lagrange_weights(x::Vector)
    n = length(x) - 1
    w = similar(x)
    w[1] = 1.0
    for j in 1:n
        for k in 0:j-1
            w[k+1] = (x[k+1] - x[j+1]) * w[k+1]
        end
        w[j+1] = prod(0:j-1) do k
            x[j+1] - x[k+1]
        end
    end
    for j in 0:n
        w[j+1] = 1/w[j+1]
    end
    return w
end

# special nodes distributions and the corresponding weights

"""
    cheb1nodes(n,a,b)

Return the `n` Chebyshev points of the first kind on the interval `[a,b]`.
"""
function cheb1nodes(n::NTuple{N},lc,uc) where {N}
    iter = cheb1nodes_iter(n,lc,uc)
    return map(x->SVector{N,Float64}(x),iter)
end
function cheb1nodes(n::NTuple{N},rec::AbstractHyperRectangle{N}) where {N}
    cheb1nodes(n,low_corner(rec),high_corner(rec))
end
function cheb1nodes(n::Integer,a::Number,b::Number)
    xcheb::Vector{Float64} = cheb1nodes(n)
    c0 = (a+b)/2
    c1 = (b-a)/2
    return @. c0 + c1*xcheb
end
function cheb1nodes(n)
    return [-cos( (2*i-1)*π / (2*n) ) for i in 1:n]
end

function cheb1nodes_iter(n::NTuple{N},lc,uc) where {N}
    nodes1d = ntuple(d->cheb1nodes(n[d],lc[d],uc[d]),N)
    return iter = Iterators.product(nodes1d...)
end
function cheb1nodes_iter(n::NTuple{N},rec::AbstractHyperRectangle{N}) where {N}
    cheb1nodes_iter(n,low_corner(rec),high_corner(rec))
end

function cheb1weights(n::NTuple{N}) where {N}
    weights1d = ntuple(d->cheb1weights(n[d]),N)
    iter = Iterators.product(weights1d...)
    return map(x->prod(x),iter)
end
function cheb1weights(i,n)
    sgn = (-1)^(i)
    val = sin((2*i-1)*π/(2*n))
    return sgn*val
end
function cheb1weights(n)
    return [(-1)^(i)*sin((2*i-1)*π/(2*n)) for i in 1:n]
end

"""
    cheb2nodes(n,a,b)

Return the `n` Chebyshev points of the second kind on the interval `[a,b]`. The
nodes are nested in the following sense: `cheb2nodes(n,a,b) ==
cheb2nodes(2n-1,a,b)[1:2:end]`.
"""
function cheb2nodes(n::NTuple{N},lc,uc) where {N}
    iter = cheb2nodes_iter(n,lc,uc)
    return map(x->SVector{N,Float64}(x),iter)
end
function cheb2nodes(n::NTuple{N},rec::AbstractHyperRectangle{N}) where {N}
    cheb2nodes(n,low_corner(rec),high_corner(rec))
end
function cheb2nodes(n::Integer,a::Number,b::Number)
    xcheb::Vector{Float64} = cheb2nodes(n)
    c0 = (a+b)/2
    c1 = (b-a)/2
    return @. c0 + c1*xcheb
end
function cheb2nodes(n)
    return [cos((i-1)*π /(n-1)) for i in 1:n]
end

function cheb2nodes_iter(n::NTuple{N},lc,uc) where {N}
    nodes1d = ntuple(d->cheb2nodes(n[d],lc[d],uc[d]),N)
    iter = Iterators.product(nodes1d...)
    return iter
end
function cheb2nodes_iter(n::NTuple{N},rec::AbstractHyperRectangle{N}) where {N}
    return cheb2nodes_iter(n,low_corner(rec),high_corner(rec))
end

function cheb2weights(n::NTuple{N}) where {N}
    weights1d = ntuple(d->cheb2weights(n[d]),N)
    iter = Iterators.product(weights1d...)
    return map(x->prod(x),iter)
end

function cheb2weights(i,n)
    sgn = (-1)^(i)
    val = i==n || i == 1 ? 0.5 : 1.0
    return sgn*val
end

function cheb2weights(n)
    return [cheb2weights(i,n) for i in 1:n]
end
