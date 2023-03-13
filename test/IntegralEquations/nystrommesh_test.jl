using Test
using LinearAlgebra
using Gmsh
import WavePropBase as WPB

@testset "Parametric surfaces" begin
    @testset "Area" begin
        @testset "Cube" begin
            # generate a mesh
            WPB.clear_entities!()
            (lx, ly, lz) = widths = (1.0, 1.0, 2.0)
            n = 2
            Ω = WPB.Domain(WPB.Box(; widths))
            Γ = WPB.boundary(Ω)
            M = WPB.meshgen(Γ, (n, n))
            msh = WPB.NystromMesh(M; qorder=1)
            A = 2 * (lx * ly + lx * lz + ly * lz)
            @test A ≈ WPB.integrate(x -> 1, msh)
        end
        @testset "Sphere" begin
            WPB.clear_entities!()
            r = 0.5
            n = 4
            Ω = WPB.Domain(WPB.Ball(; radius=r))
            Γ = WPB.boundary(Ω)
            M = WPB.meshgen(Γ, (n, n))
            msh = WPB.NystromMesh(M; qorder=4)
            area = WPB.integrate(x -> 1, msh)
            @test isapprox(area, 4 * π * r^2, atol=5e-6)
        end
        @testset "Disk" begin
            WPB.clear_entities!()
            r = rx = ry = 0.5
            n = 4
            Ω = WPB.Domain(WPB.Disk(; radius=r))
            Γ = WPB.boundary(Ω)
            M = WPB.meshgen(Γ, n)
            msh = WPB.NystromMesh(M; qorder=4)
            P = 2π * r
            @test isapprox(P, WPB.integrate(x -> 1, msh); atol=1e-6)
        end
    end
end

@testset "Gmsh" begin
    @testset "Area/volume" begin
        @testset "Cube" begin
            # generate a mesh
            (lx, ly, lz) = widths = (1.0, 1.0, 2.0)
            Ω, M = WPB.@gmsh begin
                WPB.clear_entities!()
                gmsh.model.occ.addBox(0, 0, 0, lx, ly, lz)
                gmsh.model.occ.synchronize()
                gmsh.model.mesh.generate(3)
                Ω = WPB.gmsh_import_domain(; dim=3)
                M = WPB.gmsh_import_mesh(Ω; dim=3)
                return Ω, M
            end
            ∂Ω = WPB.boundary(Ω)
            mesh = WPB.NystromMesh(view(M, ∂Ω); qorder=1)
            A = 2 * (lx * ly + lx * lz + ly * lz)
            @test A ≈ WPB.integrate(x -> 1, mesh)
            # generate a Nystrom mesh for volume
            mesh = WPB.NystromMesh(view(M, Ω); qorder=1)
            V = prod(widths)
            # sum only weights corresponding to tetras
            @test V ≈ WPB.integrate(x -> 1, mesh)
        end
        @testset "Sphere" begin
            # create a mesh using GmshSDK package
            r = 0.5
            WPB.clear_entities!()
            Ω, M = WPB.@gmsh begin
                WPB.set_meshsize(0.1)
                gmsh.model.occ.addSphere(0, 0, 0, r)
                gmsh.model.occ.synchronize()
                gmsh.model.mesh.generate(3)
                Ω = WPB.gmsh_import_domain(; dim=3)
                M = WPB.gmsh_import_mesh(Ω; dim=3)
                return Ω, M
            end
            Γ = WPB.boundary(Ω)
            mesh = WPB.NystromMesh(view(M, Γ); qorder=4) # NystromMesh of surface Γ
            area = WPB.integrate(x -> 1, mesh)
            @test isapprox(area, 4 * π * r^2, atol=5e-2)
            mesh = WPB.NystromMesh(view(M, Ω); qorder=4) # Nystrom mesh of volume Ω
            volume = WPB.integrate(x -> 1, mesh)
            @test isapprox(volume, 4 / 3 * π * r^3, atol=1e-2)
        end
        @testset "Circle" begin
            r = rx = ry = 0.5
            WPB.clear_entities!()
            Ω, M = WPB.@gmsh begin
                WPB.clear_entities!()
                gmsh.model.occ.addDisk(0, 0, 0, rx, ry)
                gmsh.model.occ.synchronize()
                gmsh.model.mesh.generate(2)
                Ω = WPB.gmsh_import_domain(; dim=2)
                M = WPB.gmsh_import_mesh(Ω; dim=2)
                return Ω, M
            end
            Γ = WPB.boundary(Ω)
            mesh = WPB.NystromMesh(view(M, Ω); qorder=2)
            A = π * r^2
            # test area
            @test isapprox(A, WPB.integrate(x -> 1, mesh); atol=1e-2)
            # test perimeter
            mesh = WPB.NystromMesh(view(M, Γ); qorder=2)
            P = 2π * r
            @test isapprox(P, WPB.integrate(x -> 1, mesh); atol=1e-2)
        end
    end
end

@testset "Near interaction list" begin
    @test_broken true == false # missing tests
end
