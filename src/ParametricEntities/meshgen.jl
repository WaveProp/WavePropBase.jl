# Basic utilities for meshing parametric entities and generate a `GenericMesh`

function meshgen!(mesh::GenericMesh,ent::ParametricEntity,sz)
    @assert ambient_dimension(ent) == ambient_dimension(mesh)
    N    = geometric_dimension(ent)
    @assert length(sz) == N
    # extract relevant fields and mesh the entity
    f    = parametrization(ent)
    d    = domain(ent)
    els  = _meshgen(f,d,sz)
    # push related information to mesh
    E      = eltype(els)
    vals   = get!(mesh.elements,E,Vector{E}())
    istart = length(vals) + 1
    append!(vals,els)
    iend   = length(vals)
    haskey(mesh.ent2tags,ent) && @warn "entity $(key(ent)) already present in mesh"
    mesh.ent2tags[ent] = Dict(E=>collect(istart:iend)) # add key
    return mesh
end

"""
    _meshgen(f,d,sz)

Create a UniformCartesianMesh` of `d` push-forward map. The cartesian mesh has
size `sz`, and is uniform in parameter space.

"""
function _meshgen(f,d,sz)
    grid = UniformCartesianMesh(d,sz)
    iter = ElementIterator(grid)
    els  = [ParametricElement(f,i) for i in iter]
    return els
end

"""
    meshgen!(mesh,Ω,sz)

Similar to [`meshgen`](@ref), but append entries to `mesh`.
"""
function meshgen!(mesh::GenericMesh,Ω::Domain,num_elements)
    @assert length(entities(Ω)) == length(num_elements)
    for (ent,sz) in zip(Ω,num_elements)
        ent isa ParametricEntity || error("meshgen! only works on parametric entites")
        meshgen!(mesh,ent,sz)
    end
    return mesh
end

"""
    meshgen(Ω::Domain,num_elements)
    meshgen(Ω::Domain;meshsize)

Generate a `GenericMesh` for the domain `Ω` with `num_elements` per entity. To
specify a different number of elements per entity, `num_elements` should be a
vector with as many elements as there are entities in `Ω`. Alternatively, a
`meshsize` can be passed.

Requires the entities forming `Ω` to be `ParametricEntity`.

!!! warning
    The quality of the generated mesh created usign `meshgen` depends on the
    quality of the underlying parametrization. For complex surfaces, you are
    better off using a proper mesher such as `gmsh`.
"""
function meshgen(Ω::Domain,num_elements::Vector)
    # extract the ambient dimension for these entities (i.e. are we in 2d or
    # 3d). Only makes sense if all entities have the same ambient dimension.
    N  = ambient_dimension(first(Ω))
    @assert all(p->ambient_dimension(p)==N,entities(Ω))
    mesh = GenericMesh{N,Float64}()
    meshgen!(mesh,Ω,num_elements) # fill in
    return mesh
end
function meshgen(Ω::Domain,num_elements::Union{Int,Tuple{Int},Tuple{Int,Int}})
    n = length(entities(Ω))
    meshgen(Ω,[num_elements for _ in 1:n])
end

# helper function to compute numbef of elements given a desired meshsize
function meshgen(Ω::Domain;meshsize)
    num_elements = []
    for ent in Ω
        l = _length_per_dimension(ent)
        n = map(i->ceil(Int,i/meshsize),l)
        push!(num_elements,n)
    end
    meshgen(Ω,num_elements)
end

# typical length of the ent across each dimension, obtained by measuring the
# length of the line passing through the center of the reference element in each
# dimension
function _length_per_dimension(ent)
    N = geometric_dimension(ent)
    d = domain(ent)
    xl,xu = low_corner(d), high_corner(d)
    ntuple(N) do dim
        f = (s) -> begin
            x = svector(i->i==dim ? s[1] : 0.5,N)
            ent(x)
        end
        integrate_fejer1d(s->integration_measure(f,s),xl[dim],xu[dim])
    end
end

function integration_measure(f,x)
    jac = jacobian(f,x)
    _integration_measure(jac)
end
function _integration_measure(jac::AbstractMatrix)
    M,N = size(jac)
    if M == N
        det(jac) |> abs # cheaper when `M=N`
    else
        transpose(jac)*jac |> det |> sqrt
    end
end

# implement the ElementIterator interface
function Base.size(iter::ElementIterator{<:ParametricElement,<:GenericMesh})
    E = eltype(iter)
    M = mesh(iter)
    els::Vector{E} = M.elements[E]
    return (length(els),)
end

function Base.getindex(iter::ElementIterator{<:ParametricElement,<:GenericMesh},i::Int)
    E = eltype(iter)
    M = mesh(iter)
    els::Vector{E} = M.elements[E]
    return els[i]
end

function Base.iterate(iter::ElementIterator{<:ParametricElement,<:GenericMesh}, state=1)
    state > length(iter) && (return nothing)
    iter[state], state + 1
end
