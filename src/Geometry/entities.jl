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
`abs(tag(Ω1))=abs(tag(Ω2))`. For `ElementaryEntity` of co-dimension one, the
sign of `tag(Ω)` is used to determine the orientation of the normal vector.

Notice that this implies `dim` and `tag` of an elementary entity should uniquely
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
    jac::SMatrix = jacobian(ent, u)
    return _normal(jac)
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
        n = SVector(t[2], -t[1])
        return n / norm(n)
    elseif M == 2 # a surface in 3d
        t₁ = jac[:, 1]
        t₂ = jac[:, 2]
        n = cross(t₁, t₂)
        return n / norm(n)
    else
        notimplemented()
    end
end

"""
    struct ElementaryEntity <: AbstractEntity

The most basic representation of an [`AbstractEntity`](@ref).

# Fields:
- `dim::UInt8`: the geometrical dimension of the entity (e.g. line has `dim=1`,
  surface has `dim=2`, etc)
- `tag::Int64`: an integer tag associated to the entity
- `boundary::Vector{AbstractEntity}`: the entities of dimension `dim-1`
  forming the entity's boundary
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

function ElementaryEntity(; boundary, dim::Int=_compute_dim_from_boundary(boundary))
    t = new_tag(dim)
    return ElementaryEntity(UInt8(dim), t, boundary)
end

function _compute_dim_from_boundary(boundary)
    dmin, dmax = extrema(geometric_dimension, boundary)
    @assert dmin == dmax "all entities in `boundary` must have the same dimension"
    return dim = dmin + 1
end

function flip_normal(e::ElementaryEntity)
    msg = "flip_normal only works for entities of co-dimension one."
    @assert ambient_dimension(ent) == geometric_dimension(ent) + 1 msg
    d = geometric_dimension(e)
    t = tag(e)
    bnd = boundary(e)
    return ElementaryEntity(d, -t, bnd)
end

function normal(ent::ElementaryEntity, u)
    s = sign(tag(ent))
    jac::SMatrix = jacobian(ent, u)
    return s * _normal(jac)
end

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
[`ElementaryEntity`](@ref) for an example.
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
