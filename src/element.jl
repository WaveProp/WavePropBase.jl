"""
    abstract type AbstractElement{D,T}

Abstract element given by an interpolation scheme of [`return_type`](@ref) `T` over
the domain `D<:ReferenceShape`.

Instances `el` of `AbstractElement` are expected to implement:
- `el(x̂)`: evaluate the interpolation scheme at the coordinate `x̂ ∈ D`.
- `jacobian(el,x̂)` : evaluate the jacobian matrix of the interpolation at the
    coordinate `x ∈ D`.

!!! note
    For performance reasons, both `el(x̂)` and `jacobian(el,x̂)` should
    typically operate on static arrays/vectors.
"""
abstract type AbstractElement{D,T} end

function (el::AbstractElement)(x)
    abstractmethod(typeof(el))
end

function jacobian(el::AbstractElement, x)
    abstractmethod(typeof(el))
end

domain(::SType{<:AbstractElement{D}}) where {D <: ReferenceShape} = D()

return_type(el::AbstractElement{D,T}) where {D,T} = T

domain_dimension(t::SType{<:AbstractElement}) = dimension(domain(t))

range_dimension(el::AbstractElement{R,T})        where {R,T} = length(T)

range_dimension(t::Type{<:AbstractElement{R,T}}) where {R,T} = length(T)

"""
    struct LagrangeElement{D,Np,T} <: AbstractElement{D,T}

# Fields:
- `vals::SVector{Np,T}`
"""
struct LagrangeElement{D,Np,T} <: AbstractElement{D,T}
    vals::SVector{Np,T}
end

vals(el::LagrangeElement) = el.vals

# a contructor which infers the extra information from nodes
function LagrangeElement{R}(vals::SVector{Np,T}) where {R,Np,T}
    LagrangeElement{R,Np,T}(vals)
end

"""
    reference_nodes(::LagrangeElement{D,Np,T})

Return the reference nodes on `D` used for the polynomial interpolation. The
function values on these nodes completely determines the interpolating
polynomial.

We use the same convention as `gmsh` for defining the reference nodes and their
order (see [node
ordering](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering) on `gmsh`
documentation).
"""
function reference_nodes end

# constructor which converts each entry to a Point, and then creates an SVector
# of that.
function LagrangeElement{R}(vals) where {R}
    vals = SVector.(vals)
    LagrangeElement{R}(SVector(vals))
end

# a convenience constructor to allow things like LagrangeLine(a,b) instead of LagrangeLine((a,b))
function LagrangeElement{R}(nodes...) where {R}
    LagrangeElement{R}(nodes)
end

degree(el::LagrangeElement) = degree(typeof(el))

"""
    degree(el::LagrangeElement)

The degree of the underlying polynomial used to represent this type of element.
"""
degree(el::Type{LagrangeElement{ReferenceLine,Np,T}}) where {Np,T} = Np - 1

function degree(::Type{LagrangeElement{ReferenceTriangle,Np,T}}) where {Np,T}
   p   = (-3 + sqrt(1 + 8 * Np)) / 2
   msg =
   "unable to determine degree for `LagrangeElement` over `ReferenceTriangle` containing Np=$(Np)  interpolation points. Need `Np = (p+1)*(p+2)/2` with `p` an integer."
   @assert isinteger(p) msg
   return Int(p)
end

function (el::LagrangeElement{ReferenceLine,2})(u)
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1]
end

function (el::LagrangeElement{ReferenceLine,3})(u)
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()
    v = vals(el)
    v[1] + (4 * v[3] - 3 * v[1] - v[2]) * u[1]  + 2 * (v[2] + v[1] - 2 * v[3]) * u[1]^2
end

function jacobian(el::LagrangeElement{ReferenceLine,3}, u)
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()
    v = vals(el)
    hcat((4 * v[3] - 3 * v[1] - v[2] + 4 * (v[2] + v[1] - 2 * v[3]) * u[1]))
end

function (el::LagrangeElement{ReferenceTriangle,3})(u)
    @assert length(u) == 2
    @assert u ∈ ReferenceTriangle()
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[3] - v[1]) * u[2]
end

function (el::LagrangeElement{ReferenceSquare,4})(u)
    @assert length(u) == 2
    @assert u ∈ ReferenceSquare()
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[4] - v[1]) * u[2] + (v[3] + v[1] - v[2] - v[4]) * u[1] * u[2]
end

function (el::LagrangeElement{ReferenceTetrahedron,4})(u)
    @assert length(u) == 3
    @assert u ∈ ReferenceTetrahedron()
    v = vals(el)
    v[1] + (v[2] - v[1]) * u[1] + (v[3] - v[1]) * u[2] + (v[4] - v[1]) * u[3]
end

function jacobian(el::LagrangeElement{ReferenceLine,2}, u)
    N = ambient_dimension(el)
    @assert length(u) == 1
    @assert u ∈ ReferenceLine()
    v = vals(el)
    return SMatrix{N,1}(v[2] - v[1])
end

function jacobian(el::LagrangeElement{ReferenceTriangle,3}, u)
    @assert length(u) == 2
    @assert u ∈ ReferenceTriangle()
    v = vals(el)
    hcat(
        (v[2] - v[1]),
        (v[3] - v[1])
    )
end

function jacobian(el::LagrangeElement{ReferenceSquare,4}, u)
    @assert length(u) == 2
    @assert u ∈ ReferenceSquare()
    v = vals(el)
    hcat(
        ((v[2] - v[1]) + (v[3] + v[1] - v[2] - v[4]) * u[2]),
        ((v[4] - v[1]) + (v[3] + v[1] - v[2] - v[4]) * u[1])
    )
end

function jacobian(el::LagrangeElement{ReferenceTetrahedron,4}, u)
    @assert length(u) == 3
    @assert u ∈ ReferenceTriangle()
    v = vals(el)
    hcat(
        (v[2] - v[1]),
        (v[3] - v[1]),
        (v[4] - v[1])
    )
end

# define some aliases used for representing elements of geometric nature.
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
    const LagrangeRectangle = LagrangeElement{ReferenceSquare}
"""
const LagrangeRectangle   = LagrangeElement{ReferenceSquare}
