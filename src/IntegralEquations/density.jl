"""
    struct NystromDensity{T} <: AbstractVector{T}

Discrete density with values `vals` on the [`qnodes`](@ref) of a [`NystromMesh`](@ref)
"""
struct NystromDensity{V,S<:NystromMesh} <: AbstractVector{V}
    vals::Vector{V}
    mesh::S
end

# AbstractArray interface
Base.size(σ::NystromDensity) = size(vals(σ))
Base.getindex(σ::NystromDensity, args...) = getindex(vals(σ), args...)
Base.setindex!(σ::NystromDensity, args...) = setindex!(vals(σ), args...)
Base.similar(σ::NystromDensity) = NystromDensity(similar(vals(σ)), mesh(σ))

vals(σ::NystromDensity) = σ.vals
mesh(σ::NystromDensity) = σ.mesh

"""
    NystromDensity(f::Function, X::NystromMesh)

Return a `NystromDensity` with values `f(q)` at the [`qnodes`](@ref) of `X`.

Note that the argument passsed to `f` is a [`QuadratureNode`](@ref), so that `f`
may depend on quantities other than the [`coords`](@ref) of the quadrature node
(such as the [`normal`](@ref) vector).

See also: [`QuadratureNode`](@ref)
"""
function NystromDensity(f::Function, X::NystromMesh)
    vals = [f(dof) for dof in qnodes(X)]
    return NystromDensity(vals, X)
end
