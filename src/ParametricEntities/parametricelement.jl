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
ambient_dimension(p::ParametricElement)   = length(return_type(p))

# constructor which infers the return type of f. Warn if inferred type is not
# bitstype
function ParametricElement(f,d)
    x = center(d)
    T = Base.promote_op(f,typeof(x))
    isbitstype(T) || (@warn "non bitstype detected for ParametricElement")
    D = domain(d) |> typeof
    return ParametricElement{D,T}((x)->f(d(x)))
end

function (el::ParametricElement)(u)
    @assert u ∈ domain(el)
    f = parametrization(el)
    return f(u)
end

# Plot recipes
@recipe function f(el::ParametricElement;npts=2)
    D = domain(el)
    grid   --> false
    aspect_ratio --> :equal
    label --> ""
    if D isa ReferenceLine
        s       = LinRange(0,1,npts)
        pts     = [el(v) for v in s]
        x       = [pt[1] for pt in pts]
        y       = [pt[2] for pt in pts]
        x,y
    elseif D isa ReferenceSquare
        seriestype := :surface
        xrange = LinRange(0,1,npts)
        yrange = LinRange(0,1,npts)
        pts    = [el((x,y)) for x in xrange, y in yrange]
        x      =  [pt[1] for pt in pts]
        y      =  [pt[2] for pt in pts]
        z      =  [pt[3] for pt in pts]
        x,y,z
    else
        notimplemented()
    end
end
