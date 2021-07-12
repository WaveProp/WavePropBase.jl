"""
    struct ChebInterp{N,Td,T}

Chebyshev interpolation over an `N`-dimensional `HyperRectangle`.

The main constructor takes an `Array{N,T}` of the function `vals` at the tensor
product grid formed by the one-dimensional Chebyshev nodes (of the second kind).

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

@inline function nodes(p::ChebInterp,d,i)
    lc = low_corner(domain(p))
    hc = high_corner(domain(p))
    a  = lc[d]
    b  = hc[d]
    sz = size(vals(p))
    n  = sz[d]
    x  = -cos((i-1)*π/(n-1))
    c0 = (a+b)/2
    c1 = (b-a)/2
    return c0+c1*x
end
function nodes(p,d)
    sz = size(vals(p))
    n  = sz[d]
    [nodes(p,d,i) for i in 1:n]
end
function nodes(p)
    sz = size(vals(p))
    N  = length(sz)
    ntuple(d->nodes(p,d),N)
end

@inline function weights(p::ChebInterp,d,i)
    sz = size(vals(p))
    n  = sz[d]
    sgn = (-1)^(i-1)
    δ   = (i==1 || i==n) ? 1/2 : 1
    return sgn*δ
end


# nodes(p::ChebInterp)   = p.nodes
# weights(p::ChebInterp) = p.weights

return_type(::ChebInterp{_,__,T}) where {_,__,T} = T

ambient_dimension(::ChebInterp{N}) where {N} = N

function (p::ChebInterp{N,Td,T})(x::SVector) where {N,Td,T}
    num = zero(T)
    den = zero(Td)
    for I in CartesianIndices(p.vals)
        wi     = one(Td)
        x_m_xi = one(Td)
        for d in 1:N
            i = I[d]
            wi     *= weights(p,d,i)
            xi = nodes(p,d,i)
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
