using Test
using LinearAlgebra
import WavePropBase as WPB
using Random
Random.seed!(1)

rtol = 5e-2

# for t in (:interior,:exterior)
t = :interior
σ = t == :interior ? 1 / 2 : -1 / 2
N = 2
pde = WPB.Laplace(;dim=N)
WPB.clear_entities!()
ent = N == 2 ? WPB.Disk() : WPB.Ball()
Ω = WPB.Domain(ent)
Γ = WPB.boundary(Ω)
meshsize = 0.01
qorder = 3
M = WPB.meshgen(Γ; meshsize)
mesh = WPB.NystromMesh(view(M, Γ); qorder)
xs = t == :interior ? ntuple(i -> 3, N) : ntuple(i -> 0.1, N)
T = WPB.default_density_eltype(pde)
c = rand(T)
u = (qnode) -> WPB.SingleLayerKernel(pde)(qnode, xs) * c
dudn = (qnode) -> WPB.AdjointDoubleLayerKernel(pde)(qnode, xs) * c
γ₀u = WPB.Density(u, mesh)
γ₁u = WPB.Density(dudn, mesh)
γ₀u_norm = norm(norm.(γ₀u, Inf), Inf)
γ₁u_norm = norm(norm.(γ₁u, Inf), Inf)
# single and double layer
S = WPB.SingleLayerOperator(pde, mesh)
Smat = Matrix(S)
D = WPB.DoubleLayerOperator(pde, mesh)
Dmat = Matrix(D)
e0 = norm(Smat * γ₁u - Dmat * γ₀u - σ * γ₀u, Inf) / γ₀u_norm
δS, δD = WPB.dim_correction(pde, mesh, mesh, Smat, Dmat)
Sdim, Ddim = Smat + δS, Dmat + δD
e1 = norm(Sdim * γ₁u - Ddim * γ₀u - σ * γ₀u, Inf) / γ₀u_norm
@show e1

qreg = WPB.etype2qrule(mesh,first(keys(mesh)))
qsing = WPB.SingularQuadratureRule(WPB.GaussLegendre(5*qorder),WPB.Kress(;order=3))

qsing = WPB.SingularQuadratureRule(WPB.G7K15(;atol=1e-6),WPB.Kress(;order=1))

qsing_dict = Dict(E=>qsing for E in keys(mesh))

# qsing = WPB.GaussLegendre(10)
tol = 3*meshsize
δS = WPB.singularquadrule_correction(S, qsing_dict; tol)
δD = WPB.singularquadrule_correction(D, qsing_dict; tol)

Sa, Da = Smat + δS, Dmat + δD
e2 = norm(Sa * γ₁u - Da * γ₀u - σ * γ₀u, Inf) / γ₀u_norm
@show e2
