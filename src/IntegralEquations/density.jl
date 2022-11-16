"""
    struct Density{T,S} <: AbstractVector{T}

Discrete density with values `vals` on an `AbstractMesh`
"""
struct Density{V,S<:AbstractMesh} <: AbstractVector{V}
    vals::Vector{V}
    mesh::S
end

# AbstractArray interface
Base.size(σ::Density) = size(vals(σ))
Base.getindex(σ::Density, args...) = getindex(vals(σ), args...)
Base.setindex!(σ::Density, args...) = setindex!(vals(σ), args...)
Base.similar(σ::Density) = Density(similar(vals(σ)), mesh(σ))

vals(σ::Density) = σ.vals
mesh(σ::Density) = σ.mesh

function Density(f::Function, X::NystromMesh)
    vals = [f(dof) for dof in qnodes(X)]
    return Density(vals, X)
end

function γ₀(f, X::NystromMesh)
    return Density(x -> f(coords(x)), X)
end

function γ₁(f, X::NystromMesh)
    return Density(dof -> f(coords(dof), normal(dof)), X)
end
