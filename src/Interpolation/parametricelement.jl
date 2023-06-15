"""
    ParametricElement{F,D,T} <: AbstractElement{D,T}

An element represented through a explicit function `f` mapping `D` into the
element.

See also: [`AbstractElement`](@ref)
"""
struct ParametricElement{F,D,T} <: AbstractElement{D,T}
    parametrization::F
end

parametrization(el::ParametricElement) = el.parametrization

geometric_dimension(p::ParametricElement) = geometric_dimension(domain(p))
ambient_dimension(p::ParametricElement) = length(return_type(p))

function geometric_dimension(::Type{ParametricElement{F,D,T}}) where {F,D,T}
    return geometric_dimension(D())
end

# constructor which infers the return type of f. Warn if inferred type is not
# bitstype
function ParametricElement(f, d::CartesianElement{N}) where {N}
    x = center(d)
    T = Base.promote_op(f, typeof(x))
    isbitstype(T) || (@warn "non bitstype detected for ParametricElement")
    D = ReferenceHyperCube{N}
    g = (x) -> f(d(x))
    return ParametricElement{typeof(g),D,T}(g)
end

function (el::ParametricElement)(u)
    f = parametrization(el)
    return f(u)
end

# by default use ForwardDiff for jacobian
function jacobian(el::ParametricElement, s::SVector)
    f = parametrization(el)
    return ForwardDiff.jacobian(f, s)
end
jacobian(el::ParametricElement, s) = jacobian(el, SVector(s))

function curvature(el::ParametricElement,s::SVector)
    @assert geometric_dimension(el) == 1
    @assert ambient_dimension(el)   == 2
    f1 = (s) -> parametrization(el)(s)[1]
    f2 = (s) -> parametrization(el)(s)[2]
    xp  = ForwardDiff.derivative(f1, s[1])
    yp  = ForwardDiff.derivative(f2, s[1])
    xpp = ForwardDiff.derivative(u->ForwardDiff.derivative(f1,u),s[1])
    ypp = ForwardDiff.derivative(u->ForwardDiff.derivative(f2,u),s[1])
    return (xp*ypp - yp*xpp) / (xp^2 + yp^2)^(3/2)
end

# convert to a LagrangeElement by computing the linear approximation of the
# parametric element
function LagrangeElement(el::ParametricElement{F,D,T}) where {F,D,T}
    f    = parametrization(el)
    d    = D()
    vals = map(f, vertices(d))
    return LagrangeElement{D}(vals)
end
