using Test
using LinearAlgebra
import WavePropBase as WPB
using Random
Random.seed!(1)

rtol = 5e-2

# for t in (:interior,:exterior)
for t in (:interior, :exterior)
    σ = t == :interior ? 1 / 2 : -1 / 2
    for N in (2, ) # 3d not implemented yet
        ops = (WPB.Laplace(; dim=N), WPB.Helmholtz(; k=1.2, dim=N))
        for pde in ops
            @testset "Greens identity ($t) $(N)d $pde" begin
                WPB.clear_entities!()
                ent = N == 2 ? WPB.Disk() : WPB.Ball()
                Ω = WPB.Domain(ent)
                Γ = WPB.boundary(Ω)
                meshsize = 0.2
                M = WPB.meshgen(Γ; meshsize)
                mesh = WPB.NystromMesh(view(M, Γ); qorder=3)
                xs = t == :interior ? ntuple(i -> 3, N) : ntuple(i -> 0.1, N)
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
                δS = WPB.hcubature_correction(S; max_dist=5*meshsize, maxevals=20)
                δD = WPB.hcubature_correction(D; max_dist=5*meshsize, maxevals=20)
                Smat, Dmat = S0 + δS, D0 + δD
                e1 = norm(Smat * γ₁u - Dmat * γ₀u - σ * γ₀u, Inf) / γ₀u_norm
                @testset "Single/double layer $(string(pde))" begin
                    @test norm(e0, Inf) > 10 * norm(e1, Inf)
                    @test norm(e1, Inf) < rtol
                end
            end
        end
    end
end
