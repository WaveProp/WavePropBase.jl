"""
    QuadratureNode{N,T<:Real}

A point in `ℝᴺ` with a `weight` for performing numerical integration.

A `QuadratureNode` can optionally store a `normal` vector and a `curvature`
scalar, given by the divergence of the normal vector.
"""
struct QuadratureNode{N,T<:Real}
    coords::SVector{N,T}
    weight::T
    normal::Union{Nothing,SVector{N,T}}
    curvature::Union{Nothing,T}
end

coords(q::QuadratureNode) = q.coords
weight(q::QuadratureNode) = q.weight
normal(q::QuadratureNode) = q.normal
curvature(q::QuadratureNode) = q.curvature
center(q::QuadratureNode) = q.coords

"""
    struct NystromMesh{N,T} <: AbstractMesh{N,T}

A mesh data structure for solving boundary integral equation using Nyström
methods.

A `NystromMesh` can be constructed from a `mesh::AbstractMesh` and a dictionary
`etype2qrule` mapping element types in `mesh` (the keys) to an appropriate
quadrature rule for integration over elements of that type (the value).

The degrees of freedom in a `NystromMesh` are associated to nodal values at the
quadrature nodes, and are represented using [`QuadratureNode`](@ref)s.
"""
Base.@kwdef struct NystromMesh{N,T} <: AbstractMesh{N,T}
    mesh::AbstractMesh{N,T} = GenericMesh{N,T}()
    etype2qrule::Dict{DataType,AbstractQuadratureRule} = Dict{DataType,AbstractQuadratureRule}()
    qnodes::Vector{QuadratureNode{N,T}} = Vector{QuadratureNode{N,T}}()
    etype2qtags::Dict{DataType,Matrix{Int}} = Dict{DataType,Matrix{Int}}()
end

# getters
mesh(m::NystromMesh) = m.mesh

"""
    qnodes(Q::NystromMesh)

Return a vector of [`QuadratureNode`](@ref)s associated to the mesh `Q`.
"""
qnodes(m::NystromMesh) = m.qnodes
etype2qrule(m::NystromMesh, E) = m.etype2qrule[E]
etype2qtags(m::NystromMesh, E::DataType) = m.etype2qtags[E]
Base.keys(m::NystromMesh) = keys(m.etype2qtags)
Base.getindex(m::NystromMesh, E::DataType) = mesh(m)[E]

Base.isempty(m::NystromMesh) = isempty(qnodes(m))

# generators for iterating over fields of dofs
qcoords(m::NystromMesh) = (coords(q) for q in qnodes(m))
qweights(m::NystromMesh) = (weight(q) for q in qnodes(m))
qnormals(m::NystromMesh) = (normal(q) for q in qnodes(m))

function Base.show(io::IO, msh::NystromMesh)
    return print(io, " NystromMesh with $(length(qnodes(msh))) quadrature nodes")
end

"""
    integrate(f,msh::NystromMesh)

Compute `∑ᵢ f(qᵢ)wᵢ`, where the `qᵢ` are the [`qnodes`](@ref)s of `msh`,
and `wᵢ` are the `qweights`.

Note that you must define `f(::QuadratureNode)`: use `coords(q)` and `normal(q)`
if you need to access the coordinate or normal vector at que quadrature node.
"""
function integrate(f, msh::NystromMesh)
    return sum(q -> f(q) * weight(q), qnodes(msh))
end

"""
    NystromMesh(msh::AbstractMesh,e2qrule::Dict)
    NystromMesh(msh::AbstractMesh;qorder)

Construct a `NystromMesh` with the quadrature `q = e2qrule[E]` applied to each
element type `E` in msh. If an `order` keyword is passed, a default quadrature
of the desired order is used for each element type.
"""
function NystromMesh(msh::AbstractMesh{N,T}, e2qrule::Dict; curvature=false) where {N,T}
    # initialize mesh with empty fields
    nys_msh = NystromMesh{N,T}(msh,
                               e2qrule,
                               QuadratureNode{N,T}[],
                               Dict{DataType,Matrix{Int}}())
    # loop element types and generate quadrature for each
    for E in keys(msh)
        iter = msh[E]
        qrule = e2qrule[E]
        # dispatch to type-stable method
        _build_nystrom_mesh!(nys_msh, iter, qrule, curvature)
    end
    return nys_msh
end

function NystromMesh(msh::AbstractMesh; qorder, curvature=false)
    e2qrule = Dict(E => qrule_for_reference_shape(domain(E), qorder) for E in keys(msh))
    return NystromMesh(msh, e2qrule; curvature)
end

@noinline function _build_nystrom_mesh!(msh, iter, qrule::AbstractQuadratureRule, curvature_)
    N = ambient_dimension(msh)
    E = eltype(iter)
    x̂, ŵ = qrule() #nodes and weights on reference element
    num_nodes = length(ŵ)
    M = geometric_dimension(domain(E))
    istart = length(qnodes(msh)) + 1
    for el in iter
        # and all qnodes for that element
        for (x̂i,ŵi) in zip(x̂,ŵ)
            x = el(x̂i)
            μ = integration_measure(el, x̂i)
            w = μ * ŵi
            ν = N - M == 1 ? normal(el,x̂i) : nothing
            κ = curvature_ ? curvature(el,x̂i) : nothing
            qnode = QuadratureNode(x, w, ν, κ)
            push!(qnodes(msh), qnode)
        end
    end
    iend = length(qnodes(msh))
    @assert !haskey(msh.etype2qtags, E)
    msh.etype2qtags[E] = reshape(collect(istart:iend), num_nodes, :)
    return msh
end

"""
    NystromMesh(Ω;meshsize,qorder)

Create a `NystromMesh` with elements of size `meshsize` and quadrature order
`qorder`.

A mesh is first generated using [`meshgen`](@ref), and then a `NystromMesh` is
created on top of it with the quadrature information.
"""
function NystromMesh(Ω::Domain;meshsize,qorder,curvature=false)
    msh = meshgen(Ω;meshsize)
    NystromMesh(msh; qorder, curvature)
end
NystromMesh(ent::AbstractEntity;kwargs...) = NystromMesh(Domain(ent);kwargs...)
