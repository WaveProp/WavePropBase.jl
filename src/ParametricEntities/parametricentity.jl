"""
    ParametricEntity <: AbstractEntity

A geometric entity given by an explicit `parametrization` field. A
`ParametricEntity` `ent` should be a callable through `ent(x)` where `x ∈
domain`.

See also: [`AbstractEntity`](@ref)
"""
struct ParametricEntity <: AbstractEntity
    dim::UInt8
    tag::Int
    parametrization::Any
    domain::Any
    boundary::Vector{<:AbstractEntity}
    function ParametricEntity(dim, tag, f, d)
        ent = new(dim, tag, f, d, AbstractEntity[])
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        global_add_entity!(ent)
        return ent
    end
end

domain(p::ParametricEntity) = p.domain
parametrization(p::ParametricEntity) = p.parametrization
geometric_dimension(p::ParametricEntity) = p.dim
tag(p::ParametricEntity) = p.tag
boundary(p::ParametricEntity) = p.boundary

function flip_normal(ent::ParametricEntity)
    @assert ambient_dimension(ent) == geometric_dimension(ent) + 1
    return ParametricEntity(geometric_dimension(ent), -tag(ent), parametrization(ent),
                            domain(ent))
end

function ParametricEntity(f, dom)
    d = geometric_dimension(dom)
    t = new_tag(d) # automatically generate a new (valid) tag
    return ParametricEntity(d, t, f, dom)
end

function return_type(p::ParametricEntity)
    # NOTE: this relies on the brittle promote_op
    d = domain(p)
    x = center(d)
    f = parametrization(p)
    T = Base.promote_op(f, typeof(x))
    return T
end

function ambient_dimension(p::ParametricEntity)
    # HACK: evaluate the parametrization at the center of the domain to get the
    # ambient dimension. This assume the parametrization is "correct" and always
    # returns a vector of the same length
    d = domain(p)
    x = center(d)
    f = parametrization(p)
    return length(f(x))
end

function (par::ParametricEntity)(x)
    # @assert x in domain(par) "x=$x ∉ $(domain(par))"
    return par.parametrization(x)
end

"""
    line(a,b)

Create a straight line connecting points `a` and `b`. Returns an instance
of [`ParametricEntity`](@ref).
"""
function line(a::SVector, b::SVector)
    f = (u) -> a + u[1] * (b - a)
    d = HyperRectangle(0.0, 1.0)
    return ParametricEntity(f, d)
end
line(a, b) = line(SVector(a), SVector(b))
