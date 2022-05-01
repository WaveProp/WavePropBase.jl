using Test
import WavePropBase as WPB
using GmshSDK
using StaticArrays

@testset "Sphere geo" begin
    # Test the simple sphere geometry
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT,"test","IO","sphere.geo")
    Ω = @gmsh begin
        gmsh.model.getEntities()
        GmshSDK.read_geo(fname)
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
    Ω, M  = @gmsh begin
         GmshSDK.read_msh(fname)
    end
end

@testset "Element iterator" begin
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT,"test","IO","sphere.msh")
    Ω, M  = @gmsh begin
         GmshSDK.read_msh(fname)
    end
    E    = keys(M) |> first
    iter = WPB.ElementIterator{E}(M)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end

@testset "Sub mesh" begin
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT,"test","IO","sphere.msh")
    Ω, M  = @gmsh begin
         GmshSDK.read_msh(fname)
    end
    subM  = WPB.SubMesh(M,WPB.external_boundary(Ω))
    E    = keys(subM) |> first
    iter = WPB.ElementIterator(subM,E)
    @test eltype(iter) == E
    @test length(iter) == size(M.elements[E],2)
end
