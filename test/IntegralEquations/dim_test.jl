using Test
using LinearAlgebra
import WavePropBase as WPB
using Random
Random.seed!(1)

rtol = 5e-2

# for t in (:interior,:exterior)
for t in (:interior, :exterior)
    σ = t == :interior ? 1 / 2 : -1 / 2
    for N in (2, 3)
        ops = (WPB.Laplace(; dim=N), WPB.Helmholtz(; k=1.2, dim=N))
        for pde in ops
            @testset "Greens identity ($t) $(N)d $pde" begin
                WPB.clear_entities!()
                ent = N == 2 ? WPB.Disk() : WPB.Ball()
                Ω = WPB.Domain(ent)
                Γ = WPB.boundary(Ω)
                M = WPB.meshgen(Γ; meshsize=0.2)
                mesh = WPB.NystromMesh(view(M, Γ); qorder=3)
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
                @testset "Single/double layer $(string(pde))" begin
                    @test norm(e0, Inf) > 10 * norm(e1, Inf)
                    @test norm(e1, Inf) < rtol
                end
                # adjoint double-layer and hypersingular
                K = WPB.AdjointDoubleLayerOperator(pde, mesh)
                Kmat = Matrix(K)
                H = WPB.HyperSingularOperator(pde, mesh)
                Hmat = Matrix(H)
                e0 = norm(Kmat * γ₁u - Hmat * γ₀u - σ * γ₁u, Inf) / γ₁u_norm
                δK, δH = WPB.dim_correction(pde, mesh, mesh, Kmat, Hmat; derivative=true)
                Kdim, Hdim = Kmat + δK, Hmat + δH
                e1 = norm(Kdim * γ₁u - Hdim * γ₀u - σ * γ₁u, Inf) / γ₁u_norm
                @testset "Adjoint double-layer/hypersingular $(string(pde))" begin
                    @test norm(e0, Inf) > 10 * norm(e1, Inf)
                    @test norm(e1, Inf) < rtol
                end
            end
        end
    end
end
