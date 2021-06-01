"""
    struct Domain

Represent a physical domain as a union of entities.

# See also: [`AbstractEntity`](@ref), [`ElementaryEntity`](@ref).
"""
struct Domain
    entities::Vector{AbstractEntity}
end
Domain(ω::AbstractEntity) = Domain([ω,])
Domain() = Domain(AbstractEntity[])

"""
    entities(Ω::Domain)

Return a vector of all elementary entities making up a domain.
"""
entities(Ω::Domain) = Ω.entities

function geometric_dimension(Ω::Domain)
    maximum(geometric_dimension(ent) for ent in entities(Ω))
end

"""
    length(Ω:::Domain)

The length of a domain corresponds to the number of elementary entities that
make it.
"""
Base.length(Ω::Domain) = length(entities(Ω))

"""
    skeleton(Ω::Domain)

Return all the boundaries of the domain, i.e. the domain's skeleton.
"""
skeleton(Ω::Domain) = union(Domain.(boundary.(entities(Ω)))...)

"""
    ===(Ω1::Domain,Ω2::Domain)

Two `Domain`s are equal if all their entities are equal (regardless of order).
"""
function Base.:(==)(Ω1::Domain, Ω2::Domain)
    return issetequal(entities(Ω1),entities(Ω2))
end

"""
    in(ω::ElementaryEntity,Ω::Domain)

Check whether an `ElementaryEntity` belongs to a `Domain` by recursively
checking whether it belongs to its boundary.
"""
function Base.in(ω::ElementaryEntity, Ω::Domain)
    in(ω, entities(Ω)) && return true
    for bnd in boundary(Ω)
        in(bnd,Ω) && (return true)
    end
    return false
end

Base.getindex(Ω::Domain, i) = entities(Ω)[i]

Base.lastindex(Ω::Domain) = length(Ω)

"""
    iterate(Ω::Domain)

Iterating over a domain means iterating over its entities.
"""
function Base.iterate(Ω::Domain, state=1)
    # Check if we are done iterating
    if state > length(Ω)
        return nothing
    end
    # Return (result, state)
    return (Ω[state], state+1)
end

Base.isempty(Ω::Domain) = length(entities(Ω)) == 0

function Base.setdiff(Ω1::Domain, Ω2::Domain)
    Domain(setdiff(entities(Ω1), entities(Ω2)))
end

function Base.union(Ω1::Domain,Ωs::Domain...)
    Domain(Vector{AbstractEntity}(unique(vcat(entities(Ω1),entities.(Ωs)...))))
end
Base.union(Ω::Domain) = Domain(unique(entities(Ω)))

"""
    assertequaldim(Ω1::Domain,Ω2::Domain)

Check that two domains have same dimension.

If one of the domain (or both) are empty, the assertion is assumed to be true.
"""
function assertequaldim(Ω1::Domain, Ω2::Domain)
    isempty(Ω1) || isempty(Ω2) && return true
    msg = "The dimension of the first domain should be equal to the dimension
    of the second domain."
    @assert geometric_dimension(Ω1) == geometric_dimension(Ω2) msg
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

function Base.push!(Ω::Domain,ent::ElementaryEntity)
    push!(entities(Ω),ent)
end

function Base.issubset(Ω1::Domain, Ω2::Domain)
    assertequaldim(Ω1, Ω2)
    return issubset(entities(Ω1), entities(Ω2))
end


"""Return the internal boundaries inside a domain."""
function internal_boundary(Ω::Domain)
    Ω1 = Domain(Ω[1])
    γ  = Domain()
    for ω2 in Ω[2:end]
        Ω2 = Domain(ω2)
        γ1 = Domain(vcat(boundary.(entities(Ω1))...))
        γ2 = Domain(vcat(boundary.(entities(Ω2))...))
        γ = union(γ, intersect(γ1, γ2))
        Ω1 = union(Ω1, Ω2)
    end
    return γ
end

"""Return the external boundaries inside a domain."""
function external_boundary(Ω::Domain)
    return setdiff(skeleton(Ω),internal_boundary(Ω))
end

"""
Return all tags of the elementary entities in the domain `Ω` corresponding to the dimension `d`.
"""
function Base.keys(Ω::Domain, d::Integer)
    if isempty(Ω)
        return Tuple{Int64,Int64}[]
    elseif d == geometric_dimension(Ω)
        return Vector{Tuple{Int64,Int64}}(vcat(key.(entities(Ω))))
    elseif d < geometric_dimension(Ω)
        return unique(keys(skeleton(Ω),d))
    else
        error("Asking for tags with dimension > dimension of domain")
    end
end
function Base.keys(Ω::Domain)
    isempty(Ω) ? Tuple{Int64, Int64}[] : keys(Ω, geometric_dimension(Ω))
end

"""
Return all tags of the elementary entities in the domain `Ω` corresponding to the dimensions contained in `dims`.
"""
function Base.keys(Ω::Domain, dims::Vector{T}) where T <: Integer
    tgs = Vector{Tuple{Int64, Int64}}(undef, 0)
    for d in dims push!(tgs, keys(Ω, d)...) end
    return tgs
end
