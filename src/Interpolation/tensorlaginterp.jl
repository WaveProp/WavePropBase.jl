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
struct TensorLagInterp{N,Td,T}
    vals::Array{T,N}
    nodes::SVector{N,Vector{Td}}
    weights::SVector{N,Vector{Td}}
end

interpolation_nodes(p::TensorLagInterp)   = p.nodes
interpolation_nodes(p::TensorLagInterp,d::Int,i::Int)   = p.nodes[d][i]
function interpolation_nodes(p::TensorLagInterp{N},I::CartesianIndex{N}) where {N}
    svector(N) do d
        interpolation_nodes(p,d,I[d])
    end
end

weights(p::TensorLagInterp) = p.weights
vals(p::TensorLagInterp)    = p.vals

return_type(::TensorLagInterp{_,__,T}) where {_,__,T} = T

ambient_dimension(::TensorLagInterp{N}) where {N} = N

Base.CartesianIndices(p::TensorLagInterp) = CartesianIndices(vals(p))

# if passed the values and nodes, construct weights automatically
function TensorLagInterp(vals::Array{T,N},nodes::SVector{N,Vector{Td}}) where {N,Td,T}
    weights = map(barycentric_lagrange_weights,nodes)
    TensorLagInterp{N,Td,T}(vals,nodes,weights)
end
TensorLagInterp(vals::Vector,nodes::Vector) = TensorLagInterp(vals,(nodes,))
TensorLagInterp(vals::Array,nodes::NTuple)  = TensorLagInterp(vals,SVector(nodes))

# special interpolation schemes may have explicit functions to compute
# interpolation nodes and weights.
function TensorLagInterp(vals::Array{T,N},rec::HyperRectangle{N},nodes_func,weights_func) where {T,N}
    lc,hc = low_corner(rec),high_corner(rec)
    sz    = size(vals)
    nodes = svector(N) do d
        a,b,n = lc[d],hc[d],sz[d]
        nodes_func[d](n,a,b)
    end
    weights = svector(N) do d
        n = sz[d]
        weights_func[d](n)
    end
    TensorLagInterp(vals,nodes,weights)
end
function TensorLagInterp(vals::Array{T,N},rec::HyperRectangle{N},nodes_func::Function,weights_func::Function) where {T,N}
    TensorLagInterp(vals,rec,ntuple(i->nodes_func,N),ntuple(i->weights_func,N))
end

# barycentric lagrange interpolation
function (p::TensorLagInterp{N,Td,T})(x::SVector) where {N,Td,T}
    num = zero(T)
    den = zero(Td)
    for I in CartesianIndices(p.vals)
        wi     = one(Td)
        x_m_xi = one(Td)
        for d in 1:N
            wi     *= p.weights[d][I[d]]
            x_m_xi *= p.nodes[d][I[d]] - x[d]
        end
        # FIXME: implicilty assumes that x is not one of the interpolation
        # nodes. Division by zero if that is the case.
        num += wi/x_m_xi*p.vals[I]
        den += wi/x_m_xi
    end
    num/den
end
(p::TensorLagInterp)(x::NTuple) = p(SVector(x))
(p::TensorLagInterp{1})(x::Number) = p((x,))

# multidimensional version of algorithm on page 504 of
# https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
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
    chebnodes(n,a,b)

Return the `n` Chebyshev poins of the first kind on the interval `[a,b]`. If `n`
is known statically, consider passing `Val(n)` instead for speed.
"""
function chebnodes(n,a,b)
    xcheb::Vector{Float64} = chebnodes(n)
    c0 = (a+b)/2
    c1 = (b-a)/2
    return @. c0 + c1*xcheb
end
function chebnodes(n)
    [-cos( (2*i-1)*π / (2*n) ) for i in 1:n]
end
@generated function chebnodes(::Val{N}) where {N}
    [-cos( (2*i-1)*π / (2*N) ) for i in 1:N]
end

function chebweight(i,n)
    sgn = (-1)^(i)
    val = sin((2*i-1)*π/(2*n))
    return sgn*val
end

function chebweights(n)
    [(-1)^(i)*sin((2*i-1)*π/(2*n)) for i in 1:n]
end

function chebweights(::Val{N}) where {N}
    [(-1)^(i)*sin((2*i-1)*π/(2*N)) for i in 1:N]
end
