"""
    abstract type AbstractEntity

Entity of geometrical nature. Identifiable throught its `(dim,tag)` key.
"""
abstract type AbstractEntity end

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

"""
    ElementaryEntity(dim,tag)

Construct an [`ElementaryEntity`](@ref) with an empty boundary .
"""
function ElementaryEntity(dim,tag)
    ElementaryEntity(dim,tag,ElementaryEntity[])
end
ElementaryEntity(dim) = ElementaryEntity(dim,_new_tag(dim))

geometric_dimension(ω::ElementaryEntity) = ω.dim

"""
    key(ω::ElementaryEntity)

Return the unique `(dim,tag)` key for the elementary entity.
"""
key(ω::ElementaryEntity) = (geometric_dimension(ω), ω.tag)

"""
    boundary(ω::ElementaryEntity)

Return the vector of elementary entities making the boundary.
"""
boundary(ω::ElementaryEntity) = ω.boundary

"""
    ==(Ω1::AbstractEntity,Ω2::AbstractEntity)

Two elementary entities are considered equal if their `dim` and `tag` fields
match.

Notice that this implies `dim` and `tag` of an elementary entity should uniquely
define it, and therefore global variables like [`TAGS`](@ref) are needed to make
sure newly created `AbstractEntities` have a new `(dim,tag)` identifier.
"""
function Base.:(==)(Ω1::AbstractEntity, Ω2::AbstractEntity)
    d1,t1 = key(Ω1)
    d2,t2 = key(Ω2)
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
    d,t = key(ent)
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

clear_tags!() = empty!(TAGS)
clear_entities!() = empty!(ENTITIES)
function clear!()
    clear_tags!()
    clear_entities!()
    nothing
end