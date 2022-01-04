"""
    abstract type AbstractElement{D,T}

Elements given by a fixed interpolation schemes mapping points on the the domain
`D<:AbstractReferenceShape` (of singleton type) to a value of type `T`.

Instances `el` of `AbstractElement` are expected to implement:
- `el(x̂)`: evaluate the interpolation scheme at the (reference) coordinate `x̂
  ∈ D`.
- `jacobian(el,x̂)` : evaluate the jacobian matrix of the interpolation at the
  (reference) coordinate `x ∈ D`.

!!! note
    For performance reasons, both `el(x̂)` and `jacobian(el,x̂)` should
    take as input a `StaticVector` and output a static vector or static array.
"""
abstract type AbstractElement{D,T} end

function (el::AbstractElement)(x)
    abstractmethod(el)
end

function jacobian(el,x)
    abstractmethod(el)
end

function WavePropBase.normal(el::AbstractElement, u)
    @assert u ∈ domain(el)
    jac = jacobian(el, u)
    normal(jac)
end

domain(::SType{<:AbstractElement{D}}) where {D <: AbstractReferenceShape} = D()

return_type(el::AbstractElement{D,T}) where {D,T} = T

domain_dimension(t::SType{<:AbstractElement}) = dimension(domain(t))

range_dimension(el::AbstractElement{R,T})        where {R,T} = length(T)

range_dimension(t::Type{<:AbstractElement{R,T}}) where {R,T} = length(T)

"""
    struct LagrangeElement{D,Np,T} <: AbstractElement{D,T}

Standard element over `D <: AbstractReferenceShape` commonly used in finite
element methods. The underlying polynomial space is [`Pk{D,K}`](@ref), and its
interpolant maps the `Np` `reference_nodes` in `D` to `Np` values of type `T`
stored in the field `vals`.
"""
struct LagrangeElement{D,Np,T} <: AbstractElement{D,T}
    vals::SVector{Np,T}
end

vals(el::LagrangeElement) = el.vals

"""
    reference_nodes(el::LagrangeElement)

Return the reference nodes on `domain(el)` used for the polynomial
interpolation. The function values on these nodes completely determines the
interpolating polynomial.

We use the same convention as `gmsh` for defining the reference nodes and their
order (see [node
ordering](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering) on `gmsh`
documentation).
"""
function reference_nodes end

# a contructor which infers the extra information from nodes
function LagrangeElement{R,K}(vals::SVector{Np,T}) where {R,K,T,Np}
    LagrangeElement{R,K,T}(vals)
end

# a contructor which infers the extra information from nodes
function LagrangeElement{D}(nodes::SVector{Np,T}) where {D,Np,T}
    LagrangeElement{D,Np,T}(nodes)
end

# constructor which converts each entry to a Point, and then creates an SVector
# of that.
function LagrangeElement{D}(vals) where {D}
    vals = SVector.(vals)
    LagrangeElement{D}(SVector(vals))
end

# a convenience constructor to allow things like LagrangeLine(a,b) instead of LagrangeLine((a,b))
function LagrangeElement{D}(nodes...) where {D}
    LagrangeElement{D}(nodes)
end

"""
    degree(el::LagrangeElement)

The polynomial degree of the element. A `LagrangeElement` of degree `K` and
domain `D` belongs to the space [`Pk{D,K}`](@ref).
"""
function degree(::SType{LagrangeElement{D,Np}})::Int where {D,Np}
    if D == ReferenceLine
        return Np-1
    elseif D == ReferenceTriangle
        K   = (-3 + sqrt(1+8*Np))/2
        return K
    elseif D == ReferenceTetrahedron
        notimplemented()
    elseif D == ReferenceSquare
        return sqrt(Np) - 1
    elseif D == ReferenceCube
        return Np^(1/3) - 1
    else
        notimplemented()
    end
end

function polynomial_space(::SType{LagrangeElement{D,Np}}) where {D,Np}
    K = degree(LagrangeElement{D,Np})
    Pk{D,K}()
end

"""
    const LagrangePoint{N,T} = LagrangeElement{ReferencePoint,1,SVector{N,T}}
"""
const LagrangePoint{N,T}        = LagrangeElement{ReferencePoint,1,SVector{N,T}}

"""
    const LagrangeLine = LagrangeElement{ReferenceLine}
"""
const LagrangeLine        = LagrangeElement{ReferenceLine}

"""
    const LagrangeTriangle = LagrangeElement{ReferenceTriangle}
"""
const LagrangeTriangle    = LagrangeElement{ReferenceTriangle}

"""
    const LagrangeTetrahedron = LagrangeElement{ReferenceTetrahedron}
"""
const LagrangeTetrahedron = LagrangeElement{ReferenceTetrahedron}

"""
    const LagrangeSquare = LagrangeElement{ReferenceSquare}
"""
const LagrangeSquare   = LagrangeElement{ReferenceSquare}

#=
Hardcode some basic elements.
TODO: Eventually this could/should be automated, at least for the LagrangeElements.
=#

# P1 for ReferenceLine
@fastmath function (el::LagrangeLine{2})(u)
    @assert u ∈ ReferenceLine()
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1]
end

@fastmath function jacobian(el::LagrangeLine{2}, u)
    @assert u ∈ ReferenceLine()
    v = vals(el)
    return hcat(v[2] - v[1])
end

# P2 for ReferenceLine
@fastmath function (el::LagrangeLine{3})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    v[1] + (4 * v[3] - 3 * v[1] - v[2]) * u[1]  + 2 * (v[2] + v[1] - 2 * v[3]) * u[1]^2
end

@fastmath function jacobian(el::LagrangeLine{3}, u)
    @assert u ∈ domain(el)
    v = vals(el)
    hcat(4 * v[3] - 3 * v[1] - v[2] + 4 * (v[2] + v[1] - 2 * v[3]) * u[1])
end

# P1 for ReferenceTriangle
@fastmath function (el::LagrangeTriangle{3})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[3] - v[1]) * u[2]
    # v[1]*(1-u[1]-u[2]) + v[2]*(u[1]) + v[3]*u[2]
end

@fastmath function jacobian(el::LagrangeTriangle{3}, u)
    @assert u ∈ domain(el)
    v   = vals(el)
    jac = hcat( v[2] - v[1],
                v[3] - v[1])
    return jac
end

# P2 for ReferenceTriangle
@fastmath function (el::LagrangeTriangle{6})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    return (1+u[2]*(-3+2u[2])+u[1]*(-3+2u[1]+4u[2]))*v[1] +
           u[1]*(-v[2]+u[1]*(2v[2]-4v[4])+4v[4]+u[2]*(-4v[4]+4v[5]-4v[6])) +
           u[2]*(-v[3]+u[2]*(2v[3]-4v[6])+4v[6])
end

@fastmath function jacobian(el::LagrangeTriangle{6}, u)
    @assert u ∈ domain(el)
    v = vals(el)
    hcat(
        (-3+4u[1]+4u[2])*v[1] - v[2] + u[1]*(4v[2]-8v[4]) + 4v[4] + u[2]*(-4v[4]+4v[5]-4v[6]),
        (-3+4u[1]+4u[2])*v[1] - v[3] + u[2]*(4v[3]-8v[6]) + u[1]*(-4v[4]+4v[5]-4v[6]) + 4v[6]
    )
end

# P3 for ReferenceTriangle
# source: https://www.math.uci.edu/~chenlong/iFEM/doc/html/dofP3doc.html
@fastmath function (el::LagrangeTriangle{10})(u)
    @assert u ∈ domain(el)
    λ₁ = 1-u[1]-u[2]
    λ₂ = u[1]
    λ₃ = u[2]
    ϕ₁ = 0.5*(3λ₁-1)*(3λ₁-2)*λ₁
    ϕ₂ = 0.5*(3λ₂-1)*(3λ₂-2)*λ₂
    ϕ₃ = 0.5*(3λ₃-1)*(3λ₃-2)*λ₃
    ϕ₄ = 4.5*λ₁*λ₂*(3λ₁-1) #894567
    ϕ₅ = 4.5*λ₁*λ₂*(3λ₂-1)
    ϕ₆ = 4.5*λ₃*λ₂*(3λ₂-1)
    ϕ₇ = 4.5*λ₃*λ₂*(3λ₃-1)
    ϕ₈ = 4.5*λ₁*λ₃*(3λ₃-1)
    ϕ₉ = 4.5*λ₁*λ₃*(3λ₁-1)
    ϕ₁₀ = 27*λ₁*λ₂*λ₃
    v = vals(el)
    return v[1]*ϕ₁ + v[2]*ϕ₂ + v[3]*ϕ₃ +
           v[4]*ϕ₄ + v[5]*ϕ₅ + v[6]*ϕ₆ +
           v[7]*ϕ₇ + v[8]*ϕ₈ + v[9]*ϕ₉ +
           v[10]*ϕ₁₀
end

@fastmath function jacobian(el::LagrangeTriangle{10,T}, u) where T
    @assert u ∈ domain(el)
    λ₁ = 1-u[1]-u[2]
    λ₂ = u[1]
    λ₃ = u[2]
    ∇λ₁ = SMatrix{1,2,eltype(T),2}(-1.,-1.)
    ∇λ₂ = SMatrix{1,2,eltype(T),2}(1.,0.)
    ∇λ₃ = SMatrix{1,2,eltype(T),2}(0.,1.)
    ∇ϕ₁ = (13.5*λ₁*λ₁-9λ₁+1)*∇λ₁
    ∇ϕ₂ = (13.5*λ₂*λ₂-9λ₂+1)*∇λ₂
    ∇ϕ₃ = (13.5*λ₃*λ₃-9λ₃+1)*∇λ₃
    ∇ϕ₄ = 4.5*((3*λ₁*λ₁-λ₁)*∇λ₂+λ₂*(6λ₁-1)*∇λ₁)
    ∇ϕ₅ = 4.5*((3*λ₂*λ₂-λ₂)*∇λ₁+λ₁*(6λ₂-1)*∇λ₂)
    ∇ϕ₆ = 4.5*((3*λ₂*λ₂-λ₂)*∇λ₃+λ₃*(6λ₂-1)*∇λ₂)
    ∇ϕ₇ = 4.5*((3*λ₃*λ₃-λ₃)*∇λ₂+λ₂*(6λ₃-1)*∇λ₃)
    ∇ϕ₈ = 4.5*((3*λ₃*λ₃-λ₃)*∇λ₁+λ₁*(6λ₃-1)*∇λ₃)
    ∇ϕ₉ = 4.5*((3*λ₁*λ₁-λ₁)*∇λ₃+λ₃*(6λ₁-1)*∇λ₁)
    ∇ϕ₁₀ = 27*(λ₁*λ₂*∇λ₃+λ₁*λ₃*∇λ₂+λ₃*λ₂*∇λ₁)
    v = vals(el)
    return v[1]*∇ϕ₁ + v[2]*∇ϕ₂ + v[3]*∇ϕ₃ +
           v[4]*∇ϕ₄ + v[5]*∇ϕ₅ + v[6]*∇ϕ₆ +
           v[7]*∇ϕ₇ + v[8]*∇ϕ₈ + v[9]*∇ϕ₉ +
           v[10]*∇ϕ₁₀
end

# P1 for ReferenceSquare
@fastmath function (el::LagrangeElement{ReferenceSquare,4})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[4] - v[1]) * u[2] + (v[3] + v[1] - v[2] - v[4]) * u[1] * u[2]
end

@fastmath function jacobian(el::LagrangeElement{ReferenceSquare,4}, u)
    @assert u ∈ domain(el)
    v = vals(el)
    hcat(
        ((v[2] - v[1]) + (v[3] + v[1] - v[2] - v[4]) * u[2]),
        ((v[4] - v[1]) + (v[3] + v[1] - v[2] - v[4]) * u[1])
    )
end

# P1 for ReferenceTetrahedron
@fastmath function (el::LagrangeElement{ReferenceTetrahedron,4})(u)
    @assert u ∈ domain(el)
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[3] - v[1]) * u[2] + (v[4] - v[1]) * u[3]
end

@fastmath function jacobian(el::LagrangeElement{ReferenceTetrahedron,4}, u)
    @assert u ∈ domain(el)
    v = vals(el)
    hcat(
        (v[2] - v[1]),
        (v[3] - v[1]),
        (v[4] - v[1])
    )
end
