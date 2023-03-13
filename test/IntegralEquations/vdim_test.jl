using Test
using LinearAlgebra
using StaticArrays
import WavePropBase as WPB
using Gmsh

pde = WPB.Laplace(; dim=3)
order = 1   # meshs order

Ψ_func = x -> cos(x[1]) * sin(x[2])
∇Ψ_func = x -> SVector(-sin(x[1]) * sin(x[2]), cos(x[1]) * cos(x[2]), 0.0)
f_func = x -> -2 * Ψ_func(x) # ΔΨ

####################################################
#    Create the domain Ω and its mesh Ωₕ
####################################################
porder = 0
h = 0.2

Ω, M = WPB.@gmsh begin
    WPB.clear_entities!()
    WPB.set_meshsize(h)
    WPB.set_meshorder(order)
    # Ω = GmshSDK.read_geo("test/kite.geo";dim=2)
    gmsh.model.occ.addSphere(0, 0, 0, 1)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    Ω = WPB.gmsh_import_domain(; dim=3)
    M = WPB.gmsh_import_mesh(Ω; dim=3)
    Ω, M
end
Γ = WPB.boundary(Ω)
Γₕ = view(M, Γ)
Ωₕ = view(M, Ω)

Ωₕ_quad = WPB.NystromMesh(Ωₕ; qorder=2)
Γₕ_quad = WPB.NystromMesh(Γₕ; qorder=5)

f = WPB.NystromDensity(dof -> f_func(dof.coords), Ωₕ_quad) # f on the volume
Ψ = WPB.NystromDensity(dof -> Ψ_func(dof.coords), Ωₕ_quad) #
γ₀Ψ = WPB.NystromDensity(dof -> Ψ_func(dof.coords), Γₕ_quad) #
γ₁Ψ = WPB.NystromDensity(dof -> dot(∇Ψ_func(dof.coords), dof.normal), Γₕ_quad) #

S₀ = Matrix(WPB.SingleLayerOperator(pde, Ωₕ_quad, Γₕ_quad))
D₀ = Matrix(WPB.DoubleLayerOperator(pde, Ωₕ_quad, Γₕ_quad))
δS, δD = WPB.dim_correction(pde, Ωₕ_quad, Γₕ_quad, S₀, D₀; location=:inside, tol=5 * h)
S = S₀ + δS
D = D₀ + δD

ref = (-D * γ₀Ψ + S * γ₁Ψ - Ψ)

V₀ = Matrix(WPB.SingleLayerOperator(pde, Ωₕ_quad))
approx0 = (V₀ * f)
er0 = ref - approx0
@show norm(er0, Inf)

δV = WPB.vdim_correction(pde, Ωₕ_quad, Ωₕ_quad, Γₕ_quad, S, D, V₀; location=:inside,
                         order=porder)
V = V₀ + δV

approx1 = (V * f)
er1 = ref - approx1
@show norm(er1, Inf)
