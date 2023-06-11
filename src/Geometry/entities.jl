"""
    abstract type AbstractEntity

Entity of geometrical nature. Identifiable throught its `(dim,tag)` key.
"""
abstract type AbstractEntity end

"""
    key(e::AbstractEntity)

The `(dim,tag)` pair used as a key to identify geometrical entities.
"""
key(e::AbstractEntity) = geometric_dimension(e), tag(e)

# reasonable defaults which assume the filds `tag` and `dim` and `boundary`
# fields exist. Some
# `AbstractEntities` need to override this method.
"""
    tag(::AbstractEntity)

Integer tag used to idetify geometrical entities.
"""
tag(e::AbstractEntity) = e.tag

"""
    geometric_dimension(x::AbstractEntity)
    geometric_dimension(Ω::Domain)

Number of degrees of freedom necessary to locally represent the geometrical
object. For example, lines have geometric dimension of 1 (whether in `ℝ²` or in
`ℝ³`), while surfaces have geometric dimension of 2.

When the argument is a `Domain`, return the largest geometric dimension
encoutered.
"""
geometric_dimension(e::AbstractEntity) = e.dim

boundary(e::AbstractEntity) = e.boundary

function Base.show(io::IO, ent::AbstractEntity)
    T = typeof(ent)
    d = geometric_dimension(ent)
    t = tag(ent)
    return print(io, "$T with (dim,tag)=($d,$t)")
end

"""
    ==(Ω1::AbstractEntity,Ω2::AbstractEntity)

Two entities are considered equal
`geometric_dimension(Ω1)==geometric_dimension(Ω2)` and
`abs(tag(Ω1))=abs(tag(Ω2))`. The sign of `tag(Ω)` is used to determine the
orientation of the normal vector.

Notice that this implies `dim` and `tag` of an entity should uniquely
define it (up to the sign of `tag`), and therefore global variables like
[`TAGS`](@ref) are needed to make sure newly created [`AbstractEntity`](@ref)
have a new `(dim,tag)` identifier.
"""
function Base.:(==)(Ω1::AbstractEntity, Ω2::AbstractEntity)
    d1, t1 = geometric_dimension(Ω1), tag(Ω1)
    d2, t2 = geometric_dimension(Ω2), tag(Ω2)
    d1 == d2 || (return false)
    abs(t1) == abs(t2) || (return false)
    # boundary(Ω1) == boundary(Ω2) || return false # this should not be needed
    return true
end
Base.hash(ent::AbstractEntity, h::UInt) = hash((geometric_dimension(ent), abs(tag(ent))), h)

function normal(ent::AbstractEntity, u)
    jac = jacobian(ent, u)
    s = sign(tag(ent))
    ν = _normal(jac)
    return s*ν
end

"""
    _normal(jac::SMatrix{M,N})

Given a an `M` by `N` matrix representing the jacobian of a codimension one
object, compute the normal vector.
"""
function _normal(jac::SMatrix{N,M}) where {N,M}
    msg = "computing the normal vector requires the element to be of co-dimension one."
    @assert (N - M == 1) msg
    if M == 1 # a line in 2d
        t = jac[:, 1] # tangent vector
        n = SVector(t[2], -t[1]) |> normalize
        return n
    elseif M == 2 # a surface in 3d
        t₁ = jac[:, 1]
        t₂ = jac[:, 2]
        n = cross(t₁, t₂) |> normalize
        return n
    else
        notimplemented()
    end
end

"""
    ParametricEntity <: AbstractEntity

A geometric entity given by an explicit `parametrization` field. A
`ParametricEntity` should be a callable through `ent(x)` where `x ∈
domain(ent)`.

See also: [`AbstractEntity`](@ref)
"""
struct ParametricEntity <: AbstractEntity
    dim::UInt8
    tag::Int
    parametrization::Any
    domain::Any
    boundary::Vector{<:AbstractEntity}
    function ParametricEntity(dim, tag, f, d, bnd = AbstractEntity[])
        @assert f(center(d)) isa SVector "parametrization must return an SVector"
        ent = new(dim, tag, f, d, bnd)
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
    msg = "flip_normal only works for entities of co-dimension one."
    @assert ambient_dimension(ent) == geometric_dimension(ent) + 1 msg
    f = parametrization(ent)
    return ParametricEntity(geometric_dimension(ent), -tag(ent), f, domain(ent))
end

"""
    ParametricEntity(f, dom::HyperRectangle[, bnd])
    ParametricEntity(f, lc, uc)

Create a [`ParametricEntity`](@ref) with parametrization `f` and domain `dom`.
If passed a tuple or vector `lc` and `uc`, the domain is taken to be a
[`HyperRectangle`](@ref) with lower corner `lc` and upper corner `uc`.

The a vector of entities composing the boundary can be passed a an optional last argument.

```jldoctest
WavePropBase.clear_entities!() # hide
WavePropBase.ParametricEntity(s->(cos(s[1]),sin(s[1])),0,2π)

# output
WavePropBase.ParametricEntity with (dim,tag)=(1,1)
```
"""
function ParametricEntity(f, dom::HyperRectangle, bnd=AbstractEntity[])
    d = geometric_dimension(dom)
    t = new_tag(d) # automatically generate a new (valid) tag
    # if return type of `f` is a Tuple, cast it as an SVector for convenience.
    if return_type(f) <: Tuple
        return ParametricEntity(d, t, s->SVector(f(s)), dom, bnd)
    else
        return ParametricEntity(d, t, f, dom, bnd)
    end
end
function ParametricEntity(f,a::SVector{N},b::SVector{N}) where {N}
    dom = HyperRectangle(a,b)
    return ParametricEntity(f, dom)
end
ParametricEntity(f,a::Tuple,b::Tuple) = ParametricEntity(f, SVector(a), SVector(b))
ParametricEntity(f,a::Vector,b::Vector) = ParametricEntity(f,Tuple(a),Tuple(b))
ParametricEntity(f,a::Number,b::Number) = ParametricEntity(f,Tuple(a),Tuple(b))

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
    @assert x in domain(par) "x=$x ∉ $(domain(par))"
    return par.parametrization(x)
end

function jacobian(psurf::ParametricEntity, s::SVector)
    return ForwardDiff.jacobian(psurf.parametrization, s)
end
jacobian(psurf::ParametricEntity, s) = jacobian(psurf, SVector(s))

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

# https://www.ljll.math.upmc.fr/perronnet/transfini/transfini.html
function transfinite_rectangle(bnd::Vector{<:ParametricEntity},tol=1e-12)
    # some checks that bnd describes the four boundaries of transfinite surface
    @assert length(bnd) == 4
    @assert norm(bnd[1](1)-bnd[2](0)) < tol &&
            norm(bnd[2](1)-bnd[3](0)) < tol &&
            norm(bnd[3](1)-bnd[4](0)) < tol &&
            norm(bnd[4](1)-bnd[1](0)) < tol
            "endpoints of boundary curves do not match"
    lc = map(b -> low_corner(domain(b)),bnd)
    hc = map(b -> high_corner(domain(b)),bnd)
    # TODO: reparametrize curves on [0,1] if that is not the case
    @assert all(x->x[1]==0,lc) && all(x->x[1]==1,hc) "domain of boundary curves is not [0,1]"
    # create transfinite interpolation
    par = _transfinite_rectangle_parametrization(bnd...)
    d   = HyperRectangle((0.,0.),(1.,1.))
    return ParametricEntity(par,d,bnd)
end
transfinite_rectangle(l1,l2,l3,l4) = transfinite_rectangle([l1,l2,l3,l4])

function _transfinite_rectangle_parametrization(c1,c2,c3,c4)
    p1 = c1.parametrization
    p2 = c2.parametrization
    p3 = c3.parametrization
    p4 = c4.parametrization
    P12 = p1(0)
    P34 = p3(0)
    P14 = p1(1)
    P32 = p3(1)
    par = x -> begin
        u,v = x[1],x[2]
        (1-v)*p1(u) + v*p3(1-u) + (1-u)*p2(v) + u*p4(1-v) -(
            (1-u)*(1-v)*P12 + u*v*P34 + u*(1-v)*P14 + (1-u)*v*P32
        )
    end
    return par
end

# TODO: add transfinite interpolation on triagnles, cubes, hexas

struct GmshEntity <: AbstractEntity
    dim::UInt8
    gmshtag::Int64
    tag::Int64
    boundary::Vector{GmshEntity}
    model::String
    function GmshEntity(d::Integer, gmshtag::Integer, model, boundary=GmshEntity[])
        msg = "an elementary entities in the boundary has the wrong dimension"
        for b in boundary
            @assert geometric_dimension(b) == d - 1 msg
        end
        tag = new_tag(d)
        ent = new(d, gmshtag, tag, boundary, model)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        global_add_entity!(ent)
        return ent
    end
end

gmshtag(e::GmshEntity) = e.gmshtag
model(e::GmshEntity) = e.model
geometric_dimension(e::GmshEntity) = e.dim
tag(e::GmshEntity) = e.tag
boundary(e::GmshEntity) = e.boundary

"""
    PointEntity{N,T} <: AbstractEntity

Zero-dimension geometrical entity. As a subtype of [`AbstractEntity`],(@ref) the
`(dim,tag)` of all created point entities get added to the global `ENTITIES`.
Intended usage is to build higher dimensionsional entities, and *not* to
represent regular points such as grid points.
"""
struct PointEntity <: AbstractEntity
    tag::Int
    coords::SVector
end

coords(p::PointEntity) = p.coords
tag(p::PointEntity) = p.tag

geometric_dimension(::PointEntity) = 0
ambient_dimension(p::PointEntity) = length(coords(p))

"""
    struct ElementaryEntity <: AbstractEntity

The most basic representation of an [`AbstractEntity`](@ref) containing a `dim`,
a `tag`, and a `boundary` field.
"""
struct ElementaryEntity <: AbstractEntity
    dim::UInt8
    tag::Int64
    boundary::Vector{<:AbstractEntity}
    function ElementaryEntity(d::Integer, t::Integer, boundary::Vector{<:AbstractEntity})
        msg = "an elementary entities in the boundary has the wrong dimension"
        for b in boundary
            @assert geometric_dimension(b) == d - 1 msg
        end
        ent = new(d, t, boundary)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        global_add_entity!(ent)
        return ent
    end
end

"""
    ElementaryEntity(dim,tag)

Construct an [`ElementaryEntity`](@ref) with an empty boundary .
"""
function ElementaryEntity(dim, tag)
    return ElementaryEntity(dim, tag, ElementaryEntity[])
end

ElementaryEntity(dim) = ElementaryEntity(dim, new_tag(dim))

function flip_normal(ent::ElementaryEntity)
    msg = "flip_normal only works for entities of co-dimension one."
    @assert ambient_dimension(ent) == geometric_dimension(ent) + 1 msg
    d = geometric_dimension(ent)
    t = tag(ent)
    bnd = boundary(ent)
    return ElementaryEntity(d, -t, bnd)
end

#####################################################################

# Variables and functions to globally keep track of entities

#####################################################################

"""
    const TAGS::Dict{Int,Vector{Int}}

Global dictionary storing the used entity tags (the value) for a given dimension
(the key).
"""
const TAGS = Dict{Int,Vector{Int}}()

"""
    const ENTITIES

Global dictionary storing the used entity tags (the value) for a given dimension
(the key).
"""
const ENTITIES = Dict{Tuple{Int,Int},AbstractEntity}()

"""
    global_add_entity!(ent::AbstractEntity)

Add `ent` to the global dictionary [`ENTITIES`](@ref) and update [`TAGS`](@ref)
with its `(dim,tag)` key. This function should be called by the inner constructor
of *every* [`AbstractEntity`](@ref); see the constructor of
[`ParametricEntity`](@ref) for an example.
"""
function global_add_entity!(ent::AbstractEntity)
    d, t = geometric_dimension(ent), tag(ent)
    _add_tag!(d, t) # add this tag to global list to make sure it is not used again
    msg = "overwriting ENTITIES: value in key ($d,$t) will be replaced"
    haskey(ENTITIES, (d, t)) && (@warn msg)
    ENTITIES[(d, t)] = ent
    return d, t
end

"""
    new_tag(dim)

Generate a unique tag for an `AbstractEntity` of dimension `dim`.

The implementation consists of adding one to the maximum value of `TAGS[dim]`

# See also: [`TAGS`](@ref).

"""
function new_tag(dim::Integer)
    if !haskey(TAGS, dim)
        return 1
    else
        tnew = maximum(TAGS[dim]) + 1
        return tnew
    end
end

function _add_tag!(dim, tag)
    if is_new_tag(dim, tag)
        # now add key
        if haskey(TAGS, dim)
            push!(TAGS[dim], tag)
        else
            TAGS[dim] = [tag]
        end
    else
        # print warning but don't add duplicate tag
        msg = "entity of dimension $dim and tag $tag already exists in TAGS.
       Creating a possibly duplicate entity."
        @warn msg
    end
    return TAGS
end

function is_new_tag(dim, tag)
    if haskey(TAGS, dim)
        existing_tags = TAGS[dim]
        if in(tag, existing_tags)
            return false
        end
    end
    return true
end

"""
    clear_entities!()

Empty the global variables used to keep track of the various entities
created.

# See also: [`ENTITIES`](@ref), [`TAGS`](@ref)
"""
function clear_entities!()
    empty!(TAGS)
    empty!(ENTITIES)
    return nothing
end
