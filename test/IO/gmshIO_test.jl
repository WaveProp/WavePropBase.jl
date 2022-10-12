using Test
import WavePropBase as WPB
using Gmsh
using StaticArrays

@testset "Basic tests" begin
    WPB.@gmsh begin
        WPB.clear_entities!()
        # try some basic things just to make sure it does not error
        gmsh.model.geo.addPoint(0,0,0)    # test native CAO
        gmsh.model.occ.addSphere(0,0,0,1) # test occ is available
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.partition(5) # test metis is available
        @test true == true #
    end
    # read domain from a model
    WPB.@gmsh begin
        WPB.clear_entities!()
        gmsh.model.occ.addSphere(0,0,0,1)
        gmsh.model.occ.synchronize()
        Ω = WPB.gmsh_import_domain()
        @test WPB.geometric_dimension(Ω) == 3
    end
    # read domain and mesh from a model
    Ω,M = WPB.@gmsh begin
        WPB.clear_entities!()
        gmsh.model.occ.addDisk(0,0,0,1,1)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
        Ω = WPB.gmsh_import_domain(dim=2)
        M = WPB.gmsh_import_mesh(Ω,dim=2)
        return Ω,M
    end
    Ω,M = WPB.@gmsh begin
        WPB.clear_entities!()
        gmsh.model.occ.addRectangle(0,0,0,1,1)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
        Ω = WPB.gmsh_import_domain(dim=2)
        M = WPB.gmsh_import_mesh(Ω,dim=2)
        return Ω,M
    end
    # read domain from a file
    Ω = WPB.@gmsh begin
        WPB.clear_entities!()
        gmsh.initialize()
        fname = joinpath(WPB.PROJECT_ROOT,"test","IO","circle.geo")
        WPB.gmsh_read_geo(fname;dim=2)
    end
    # read domain and mesh from a file
    Ω, msh = WPB.@gmsh begin
        WPB.clear_entities!()
        fname = joinpath(WPB.PROJECT_ROOT,"test","IO","circle.msh")
        WPB.gmsh_read_msh(fname;dim=2)
    end
end

@testset "Sphere geo" begin
    # Test the simple sphere geometry
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT,"test","IO","sphere.geo")
    Ω = WPB.@gmsh begin
        gmsh.model.getEntities()
        WPB.gmsh_read_geo(fname)
    end
    @test length(Ω) == 1
    @test WPB.geometric_dimension(Ω) == 3

    Γi = WPB.internal_boundary(Ω)
    Γe = WPB.external_boundary(Ω)

    @test length(Γi) == 0
    @test length(Γe) == 1
    @test WPB.geometric_dimension(Γe) == 2

    Ce = WPB.external_boundary(Γe)
    Ci = WPB.internal_boundary(Γe)
    S  = WPB.skeleton(Γe)
    @test length(Ce) == 2
    @test WPB.geometric_dimension(Ce) == 1
end

@testset "Sphere msh" begin
    # Test the simple sphere geometry
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT,"test","IO","sphere.msh")
    Ω, M  = WPB.@gmsh begin
         WPB.gmsh_read_msh(fname)
    end
end

@testset "Element iterator" begin
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT,"test","IO","sphere.msh")
    Ω, M  = WPB.@gmsh begin
         WPB.gmsh_read_msh(fname)
    end
    E    = keys(M) |> first
    iter = WPB.ElementIterator{E}(M)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end

@testset "Sub mesh" begin
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT,"test","IO","sphere.msh")
    Ω, M  = WPB.@gmsh begin
         WPB.gmsh_read_msh(fname)
    end
    subM  = WPB.SubMesh(M,WPB.external_boundary(Ω))
    E    = keys(subM) |> first
    iter = WPB.ElementIterator(subM,E)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end

@testset "Mix ParametricEntities and GmshEntity" begin
    WPB.clear_entities!()
    ball = WPB.Ball(;center=(-1,0,0),radius=1) # abstract entity
    Γ    = WPB.boundary(ball) |> WPB.Domain
    M    = WPB.meshgen(Γ;meshsize=0.4)
    Ω = WPB.@gmsh begin
        gmsh.model.occ.addSphere(0,0,0,1)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.model.occ.synchronize()
        Ω = WPB.gmsh_import_domain(dim=2)
        WPB.gmsh_import_mesh!(M,Ω)
        return Ω
    end
end
