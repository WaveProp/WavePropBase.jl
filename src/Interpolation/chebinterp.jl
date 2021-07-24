"""
    struct ChebInterp{N,Td,T}

Chebyshev interpolation over an `N`-dimensional `HyperRectangle`.

The main constructor takes an `Array{N,T}` of the function `vals` at the tensor
product grid formed by the one-dimensional Chebyshev nodes (of the first kind).

# Examples:
```julia
nx = 10
ny = 12
x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
f   = (x) -> cos(x[1]*x[2])
vals = [f((x,y)) for x in x, y in y]
p   = ChebInterp(SVector(x,y),vals)
p((0.1,0.2)) ≈ f((0.1,0.2))
```
"""

struct ChebInterp{N,Td,T}
    vals::Array{T,N}
    domain::HyperRectangle{N,Td}
end

vals(p::ChebInterp)    = p.vals
domain(p::ChebInterp)  = p.domain

return_type(::ChebInterp{_,__,T}) where {_,__,T} = T

ambient_dimension(::ChebInterp{N}) where {N} = N

"""
    chebnode(i,n,a=-1,b=1)

The `i`-th Chebyshev point of the first kind for a chebyshev polynomial of
degree `n-1` on the interval [a,b].
"""
function chebnode(i,n,a=-1,b=1)
    x =  -cos( (2*i-1)*π / (2*n) )
    c0 = (a+b)/2
    c1 = (b-a)/2
    return c0+c1*x
end

function chebweight(i,n)
    sgn = (-1)^(i)
    val = sin((2*i-1)*π/(2*n))
    return sgn*val
end

@inline function weights(p::ChebInterp,d,i)
    sz = size(vals(p))
    n  = sz[d]
    return chebweight(i,n)
end

@inline function interpolation_nodes(p::ChebInterp,d,i)
    lc = low_corner(domain(p))
    hc = high_corner(domain(p))
    a  = lc[d]
    b  = hc[d]
    sz = size(vals(p))
    n  = sz[d]
    chebnode(i,n,a,b)
end
function interpolation_nodes(p::ChebInterp{N},I::CartesianIndex{N}) where {N}
    svector(N) do d
        interpolation_nodes(p,d,I[d])
    end
end
function interpolation_nodes(p,d)
    sz = size(vals(p))
    n  = sz[d]
    [interpolation_nodes(p,d,i) for i in 1:n]
end
function interpolation_nodes(p)
    sz = size(vals(p))
    N  = length(sz)
    ntuple(d->interpolation_nodes(p,d),N)
end

function (p::ChebInterp{N,Td,T})(x::SVector) where {N,Td,T}
    num = zero(T)
    den = zero(Td)
    for I in CartesianIndices(p.vals)
        wi     = one(Td)
        x_m_xi = one(Td)
        for d in 1:N
            i = I[d]
            wi *= weights(p,d,i)
            xi = interpolation_nodes(p,d,i)
            x_m_xi *= xi - x[d]
        end
        # FIXME: implicilty assumes that x is not one of the interpolation
        # nodes. Division by zero if that is the case.
        num += wi/x_m_xi*p.vals[I]
        den += wi/x_m_xi
    end
    num/den
end
(p::ChebInterp)(x::NTuple) = p(SVector(x))
(p::ChebInterp{1})(x::Number) = p((x,))


# """
#     chebnode2(i,n,a=-1,b=1)

# The `i`-th Chebyshev point of the second kind for a chebyshev polynomial of
# degree `n-1` on the interval [a,b].
# """
# function chebnode2(i,n,a=-1,b=1)
#     x =  -cos((i-1)*π/(n-1))
#     c0 = (a+b)/2
#     c1 = (b-a)/2
#     return c0+c1*x
# end

# @inline function weights(p::ChebInterp,d,i)
#     sz = size(vals(p))
#     n  = sz[d]
#     sgn = (-1)^(i-1)
#     δ   = (i==1 || i==n) ? 1/2 : 1
#     return sgn*δ
# end
