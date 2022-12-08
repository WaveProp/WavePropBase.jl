"""
    ParametricElement{F,D,T} <: AbstractElement{D,T}

An element represented through a explicit function `f` mapping `D` into the
element.

See also: [`AbstractElement`](@ref)
"""
struct ParametricElement{F,D,T} <: AbstractElement{D,T}
    parametrization::F
    function ParametricElement{D,T}(f::F) where {F,D,T}
        return new{F,D,T}(f)
    end
end

parametrization(el::ParametricElement) = el.parametrization
domain(::ParametricElement{F,D,T}) where {F,D,T} = D()
return_type(::ParametricElement{F,D,T}) where {F,D,T} = T

geometric_dimension(p::ParametricElement) = geometric_dimension(domain(p))
ambient_dimension(p::ParametricElement) = length(return_type(p))

geometric_dimension(::Type{ParametricElement{F,D,T}}) where {F,D,T} = geometric_dimension(D())

# constructor which infers the return type of f. Warn if inferred type is not
# bitstype
function ParametricElement(f, d)
    x = center(d)
    T = Base.promote_op(f, typeof(x))
    isbitstype(T) || (@warn "non bitstype detected for ParametricElement")
    D = typeof(domain(d))
    return ParametricElement{D,T}((x) -> f(d(x)))
end

function (el::ParametricElement)(u)
    f = parametrization(el)
    return f(u)
end
