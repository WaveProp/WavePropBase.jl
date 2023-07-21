"""
    struct Domain

Represents a physical domain as a union of [`AbstractEntity`](@ref) objects.
"""
struct Domain
    entities::OrderedSet{AbstractEntity}
end

Domain(args...) = Domain(OrderedSet(args)) # allow for Domain(ω₁,ω₂,...) syntax
Domain(ω::AbstractEntity) = Domain([ω])
Domain(v::Vector) = Domain(OrderedSet(v))

function ambient_dimension(Ω::Domain)
    l, u = extrema(ambient_dimension(ent) for ent in entities(Ω))
    @assert l == u "ambient dimension of entities in a domain not equal"
    return u
end

"""
    entities(Ω::Domain)

Return a vector of all elementary entities making up a domain.
"""
entities(Ω::Domain) = Ω.entities

function Base.show(io::IO, d::Domain)
    ents = entities(d)
    n = length(entities(d))
    n == 1 ? print(io, "Domain with $n entity:\n") : print(io, "Domain with $n entities:")
    for ent in ents
        print(io, "\n\t $(ent)")
    end
    return io
end

function geometric_dimension(Ω::Domain)
    l, u = extrema(geometric_dimension(ent) for ent in entities(Ω))
    @assert l == u "geometric dimension of entities in a domain not equal"
    return u
end

"""
    iterate(Ω::Domain)

Iterating over a domain means iterating over its entities.
"""
Base.iterate(Ω::Domain, state=1) = iterate(entities(Ω), state)

Base.isempty(Ω::Domain) = isempty(entities(Ω))

"""
    skeleton(Ω::Domain)

Return all the boundaries of the domain, i.e. the domain's skeleton.
"""
function skeleton(Ω::Domain)
    ents = OrderedSet{AbstractEntity}()
    for ent in entities(Ω)
        union!(ents, boundary(ent))
    end
    return Domain(ents)
end

"""
    ===(Ω1::Domain,Ω2::Domain)

Two `Domain`s are equal if all their entities are equal (regardless of order).
"""
function Base.:(==)(Ω1::Domain, Ω2::Domain)
    return issetequal(entities(Ω1), entities(Ω2))
end

"""
    in(ω::ElementaryEntity,Ω::Domain)

Check whether an `ElementaryEntity` belongs to a `Domain`. This will first check
if the entity is in the domain's list of entities. If not, it will check if the
entity is in the boundary of any of the entities in the domain.
"""
function Base.in(ω::ElementaryEntity, Ω::Domain)
    ents = entities(Ω)
    Γ = skeleton(Ω)
    if in(ω, ents) || in(ω,Γ) # recurse on boundary
        return true
    else
        return false
    end
end

function Base.setdiff(Ω1::Domain, Ω2::Domain)
    return Domain(setdiff(entities(Ω1), entities(Ω2)))
end

function Base.union(Ω1::Domain, Ω2::Domain)
    ents = union(entities(Ω1), entities(Ω2))
    return Domain(ents)
end
function Base.union!(Ω1::Domain, Ω2::Domain)
    union!(entities(Ω1), entities(Ω2))
    return Ω1
end


"""
    assertequaldim(Ω1::Domain,Ω2::Domain)

Check that two domains have same dimension.

If one of the domain (or both) are empty, the assertion is assumed to be true.
"""
function assertequaldim(Ω1::Domain, Ω2::Domain)
    if isempty(Ω1) || isempty(Ω2)
        return true
    else
        msg = "The dimension of the first domain should be equal to the dimension
        of the second domain."
        @assert geometric_dimension(Ω1) == geometric_dimension(Ω2) msg
    end
end

function Base.intersect(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    Ωinter = Domain(intersect(entities(Ω1), entities(Ω2)))
    if isempty(Ωinter)
        if isempty(Ω1) || isempty(Ω2)
            return Domain()
        else
            return intersect(skeleton(Ω1), skeleton(Ω2))
        end
    else
        return Ωinter
    end
end

function Base.push!(Ω::Domain, ent::AbstractEntity)
    return push!(entities(Ω), ent)
end

function Base.issubset(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    return issubset(entities(Ω1), entities(Ω2))
end

"""
    internal_boundary(Ω::Domain)

Return the internal boundaries of a `Domain`. These are entities in
`skeleton(Ω)` which appear at least twice as a boundary of entities in `Ω`.
"""
function internal_boundary(Ω::Domain)
    seen     = OrderedSet{AbstractEntity}()
    repeated = OrderedSet{AbstractEntity}()
    for ω in entities(Ω)
        for γ in boundary(ω)
            in(γ, seen) ? push!(repeated, γ) : push!(seen, γ)
        end
    end
    return Domain(repeated)
end

"""
    external_boundary(Ω::Domain)

Return the external boundaries inside a domain. These are entities in the
skeleton of Ω which are not in the internal boundaries of Ω.
"""
function external_boundary(Ω::Domain)
    return setdiff(skeleton(Ω), internal_boundary(Ω))
end

"""
    boundary(Ω::Domain)

Return the [`external_boundary`](@ref) of a domain.
"""
boundary(Ω::Domain) = external_boundary(Ω)

"""
    flip_normal(Γ::Domain)

Reverse the orientation of the normal vector in the entities of a domain.
"""
function flip_normal(Γ::Domain)
    return Domain(flip_normal.(entities(Γ)))
end
