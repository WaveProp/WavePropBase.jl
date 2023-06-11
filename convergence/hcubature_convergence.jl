using Test
using LinearAlgebra
import WavePropBase as WPB
using Random
using Gmsh
Random.seed!(1)

atol = 1e-10
t = :interior
σ = t == :interior ? 1 / 2 : -1 / 2
N = 2
pde = WPB.Laplace(; dim=N)
pde = WPB.Helmholtz(; dim=N, k=20)
@info "Greens identity ($t) $(N)d $pde"
WPB.clear_entities!()
center = WPB.Point2D(0.1,0.3)
radius = 0.5
Ω = WPB.Disk(;center,radius) |> WPB.Domain

# Ω = WPB.Polygon(;vertices = [(1,-d/3),(1,-2d/3), (2,-2d/3), (2,-d/3)]) |> WPB.Domain
Γ = WPB.boundary(Ω)

n = 4
qorder = 2n-1

hh  = [1/2^i for i in 1:6]
ee0 = []
ee1 = []

for h in hh
    meshsize = h
    # set meshsize in gmsh
    # gmsh.initialize()
    # gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    # gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    # gmsh.model.mesh.setOrder(2)
    # WPB.clear_entities!()
    # gmsh.model.occ.addDisk(0, 0, 0, 1, 1)
    # gmsh.model.occ.synchronize()
    # gmsh.model.mesh.generate(2)
    # Ω = WPB.gmsh_import_domain(; dim=2)
    # M = WPB.gmsh_import_mesh(Ω; dim=2)
    # gmsh.finalize()
    # mesh = WPB.NystromMesh(view(M, Γ); qorder)
    mesh = WPB.NystromMesh(Γ; qorder, meshsize)
    xs = ent.center
    xs = t == :interior ? ent.center + WPB.svector(i->2*ent.radius,N) : ent.center + WBP.svector(i->0.5*ent.radius,N)
    T = WPB.default_density_eltype(pde)
    c = rand(T)
    u = (qnode) -> WPB.SingleLayerKernel(pde)(qnode, xs) * c
    dudn = (qnode) -> WPB.AdjointDoubleLayerKernel(pde)(qnode, xs) * c
    γ₀u = WPB.NystromDensity(u, mesh)
    γ₁u = WPB.NystromDensity(dudn, mesh)
    γ₀u_norm = norm(norm.(γ₀u, Inf), Inf)
    γ₁u_norm = norm(norm.(γ₁u, Inf), Inf)
    # single and double layer
    S = WPB.SingleLayerOperator(pde, mesh)
    S0 = Matrix(S)
    D = WPB.DoubleLayerOperator(pde, mesh)
    D0 = Matrix(D)
    e0 = norm(S0 * γ₁u - D0 * γ₀u - σ * γ₀u, Inf) / γ₀u_norm
    δS = WPB.hcubature_correction(S; maxdist=5 * meshsize, atol, maxevals=2000)
    δD = WPB.hcubature_correction(D; maxdist=5 * meshsize, atol, maxevals=2000)
    Smat, Dmat = S0 + δS, D0 + δD
    e1 = norm(Smat * γ₁u - Dmat * γ₀u - σ * γ₀u, Inf) / γ₀u_norm
    push!(ee0,e0)
    push!(ee1,e1)
end

##
order = n + 1
fig = plot(hh,ee0,m=:x,label="no correction",yscale=:log10,xscale=:log10)
plot!(fig,hh,ee1,m=:x,label="hcubature correction")

ref = hh.^order
iref = length(ref)
plot!(fig,hh,ee1[iref]/ref[iref]*ref,label=L"\mathcal{O}(h^%$order)",ls=:dash)
