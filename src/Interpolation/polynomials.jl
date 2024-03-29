"""
    abstract type AbstractPolynomialSpace{D}

A polynomial space over `D`. This is a vector space under polynomial addition
and scalar multiplication.
"""
abstract type AbstractPolynomialSpace{D} end

"""
    struct PolynomialSpace{D,K} <: AbstractPolynomialSpace{D}

The space of all polynomials of degree `≤K`, commonly referred to as `ℙₖ`.

The type parameter `D`, of singleton type, is used to determine the reference
domain of the polynomial basis. In particular, when `D` is a hypercube in `d`
dimensions, the precise definition is `ℙₖ = span{𝐱ᶿ : 0≤max(θ)≤ K}`; when `D` is
a `d`-dimensional simplex, the space is `ℙₖ = span{𝐱ᶿ : 0≤sum(θ)≤ K}`, where `θ ∈
𝐍ᵈ` is a multi-index.

# See also: [`monomial_basis`](@ref), [`lagrange_basis`](@ref)
"""
struct PolynomialSpace{D,K} <: AbstractPolynomialSpace{D} end
PolynomialSpace(d::AbstractReferenceShape, k::Int) = PolynomialSpace{typeof(d),k}()

function Base.show(io::IO, pk::PolynomialSpace{D,K}) where {D,K}
    return print(io, "ℙ$K : space of all polynomials over $D of degree ≤ $K")
end

"""
    dimension(space)

The length of a basis for `space`; i.e. the number of linearly independent elements
required to span `space`.
"""
function dimension(::SType{PolynomialSpace{D,K}}) where {D,K}
    if D == ReferenceLine
        return K + 1
    elseif D == ReferenceTriangle
        return (K + 1) * (K + 2) ÷ 2
    elseif D == ReferenceTetrahedron
        return (K + 1) * (K + 2) * (K + 3) ÷ 6
    elseif D == ReferenceSquare
        return (K + 1)^2
    elseif D == ReferenceCube
        return (K + 1)^3
    else
        notimplemented()
    end
end

"""
    monomial_basis(sp::PolynomialSpace)

Return an `NTuple` containing a basis of monomials `𝐱ᶿ` spanning the polynomial
space [`PolynomialSpace`](@ref).
"""
function monomial_basis end

@generated function monomial_basis(::PolynomialSpace{ReferenceLine,K}) where {K}
    # NOTE: it would be more efficient to interpolate the variable i directly
    # here.
    # the K+1 monomials x^0, x^1, ..., x^K
    b = ntuple(K + 1) do i
        return x -> x^(i - 1)
    end
    return :($b)
end

@generated function monomial_basis(::PolynomialSpace{ReferenceSquare,K}) where {K}
    # the K+1 monomials x^(0,0), x^(0,1),x^(1,0), ..., x^(K,K)
    I = CartesianIndices((K + 1, K + 1)) .- CartesianIndex(1, 1)
    N = length(I)
    b = ntuple(N) do i
        θ = Tuple(I[i]) # map linear to cartesian index
        return x -> prod(x .^ θ)
    end
    return :($b)
end

@generated function monomial_basis(::PolynomialSpace{ReferenceTriangle,K}) where {K}
    # the (K+1)*(K+2)/2 monomials x^(a,b) with a+b ≤ K
    # construct first the indices for the square, then filter only those for
    # which the sum is less than K.
    I = CartesianIndices((K + 1, K + 1)) .- CartesianIndex(1, 1)
    I = filter(I) do idx
        return sum(Tuple(idx)) ≤ K
    end
    N = (K + 1) * (K + 2) ÷ 2
    b = ntuple(N) do i
        θ = Tuple(I[i])
        return x -> prod(x .^ θ)
    end
    return :($b)
end

"""
    lagrange_basis(nodes,[sp::AbstractPolynomialSpace])

Return the set of `n` polynomials in `sp` taking the value of `1` on node `i`
and `0` on nodes `j ≂̸ i` for `1 ≤ i ≤ n`. For `N`-dimensional
tensor-product nodes represented in the form of an `SVector{N,Vector{T}}`, the
argument `sp` may be ommited.

!!! danger
    It is assumed that the value of a function on `nodes` uniquely determine a
    polynomial in `sp`.
"""
function lagrange_basis(nodes, sp::AbstractPolynomialSpace)
    N = dimension(sp)
    @assert length(nodes) == N
    basis = monomial_basis(sp)
    # compute the matrix of coeffcients of the lagrange polynomials over the
    # monomomial basis
    V = [p(x) for x in nodes, p in basis]
    C = V \ I
    lag_basis = ntuple(N) do j
        x -> sum(1:N) do i
            return C[i, j] * basis[i](x)
        end
    end
    return lag_basis
end

function lagrange_basis(nodes::Vector{<:Number})
    l = x -> prod(nodes) do xi
        return x - xi
    end
    w = barycentric_lagrange_weights(nodes)
    map(nodes, w) do xj, wj
        return x -> x == xj ? 1.0 : l(x) * wj / (x - xj)
    end
end
# other possible one-dimenional syntaxes
lagrange_basis(nodes::SVector{1,<:Vector}) = lagrange_basis(nodes[1])
function lagrange_basis(nodes::SVector{<:Any,SVector{1,T}}) where {T}
    x = Vector(reinterpret(T, nodes))
    return lagrange_basis(x)
end
