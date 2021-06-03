"""
    abstract type AbstractEntity

Entity of geometrical nature. Identifiable throught its `(dim,tag)` key.
"""
abstract type AbstractEntity end

"""
    key(e::AbstractEntity)

The `(dim,tag)` pair used as a key to identify various abstract entities.
"""
key(e::AbstractEntity) = geometric_dimension(e),tag(e)

"""
    struct ElementaryEntity <: AbstractEntity

The most basic representation of an [`AbstractEntity`](@ref).

# Fields:
- `dim::UInt8`: the geometrical dimension of the entity (e.g. line has `dim=1`,
  surface has `dim=2`, etc)
- `tag::Int64`: an integer tag associated to the entity
- `boundary::Vector{ElementaryEntity}`: the entities of dimension `dim-1`
  forming the entity's boundary
"""
struct ElementaryEntity <: AbstractEntity
    dim::UInt8
    tag::Int64
    boundary::Vector{ElementaryEntity}
    function ElementaryEntity(d::Integer, t::Integer, boundary::Vector{ElementaryEntity})
        msg = "an elementary entities in the boundary has the wrong dimension"
        for b in boundary
            @assert geometric_dimension(b) == d-1 msg
        end
        ent = new(d, t, boundary)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        _global_add_entity!(ent)
        return ent
    end
end

geometric_dimension(ω::ElementaryEntity) = ω.dim
tag(ω::ElementaryEntity) = ω.tag
boundary(ω::ElementaryEntity) = ω.boundary

"""
    ElementaryEntity(dim,tag)

Construct an [`ElementaryEntity`](@ref) with an empty boundary .
"""
function ElementaryEntity(dim,tag)
    ElementaryEntity(dim,tag,ElementaryEntity[])
end

ElementaryEntity(dim) = ElementaryEntity(dim,_new_tag(dim))

"""
    ==(Ω1::AbstractEntity,Ω2::AbstractEntity)

Two elementary entities are considered equal
`geometric_dimension(Ω1)==geometric_dimension(Ω2)` and
`abs(tag(Ω1))=abs(tag(Ω2))`. The sign of `tag(Ω)` is used to determine its
orientation.

Notice that this implies `dim` and `tag` of an elementary entity should uniquely
define it (up to the sign of `tag`), and therefore global variables like
[`TAGS`](@ref) are needed to make sure newly created `AbstractEntities` have a
new `(dim,tag)` identifier.
"""
function Base.:(==)(Ω1::AbstractEntity, Ω2::AbstractEntity)
    d1,t1 = geometric_dimension(Ω1),tag(Ω1)
    d2,t2 = geometric_dimension(Ω1),tag(Ω1)
    d1 == d2  || return false
    abs(t1) == abs(t2) || return false
    # boundary(Ω1) == boundary(Ω2) || return false # this should not be needed
    return true
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

function _global_add_entity!(ent::AbstractEntity)
    d,t = geometric_dimension(ent), tag(ent)
    _add_tag!(d,t) # add this tag to global list to make sure it is not used again
    msg = "overwriting ENTITIES: value in key ($d,$t) will be replaced"
    haskey(ENTITIES,(d,t)) && (@warn msg)
    ENTITIES[(d,t)] = ent
    return d,t
end

"""
    _new_tag(dim)

Generate a unique tag for an `AbstractEntity` of dimension `dim`.

The implementation consists of adding one to the maximum value of `TAGS[dim]`

# See also: [`TAGS`](@ref).

"""
function _new_tag(dim)
    if !haskey(TAGS,dim)
        return 1
    else
        tnew = maximum(TAGS[dim]) + 1
        return tnew
    end
end

function _add_tag!(dim,tag)
    if is_new_tag(dim,tag)
        # now add key
        if haskey(TAGS,dim)
            push!(TAGS[dim],tag)
        else
            TAGS[dim] = [tag,]
        end
    else
        # print warning but don't add duplicate tag
        msg  = "entity of dimension $dim and tag $tag already exists in TAGS.
        Creating a possibly duplicate entity."
        @warn msg
    end
    return TAGS
end

function is_new_tag(dim,tag)
    if haskey(TAGS,dim)
        existing_tags = TAGS[dim]
        if in(tag,existing_tags)
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
    nothing
end
