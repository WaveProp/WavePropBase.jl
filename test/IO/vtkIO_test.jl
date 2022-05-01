using Test
using WriteVTK
using GmshSDK
import WavePropBase as WPB

@testset "VTK export" begin
    # This test should simply not throw an error
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT,"test","IO","sphere.msh")
    Ω, M  = @gmsh begin
         GmshSDK.read_msh(fname)
    end
    WPB.vtk_mesh_file(M, joinpath(WPB.PROJECT_ROOT,"test","IO","ball")) |> vtk_save
    rm(joinpath(WPB.PROJECT_ROOT,"test","IO","ball.vtu"))
    WPB.vtk_mesh_file(M, Ω, joinpath(WPB.PROJECT_ROOT,"test","IO","ball")) |> vtk_save
    rm(joinpath(WPB.PROJECT_ROOT,"test","IO","ball.vtu"))
    WPB.vtk_mesh_file(M, WPB.external_boundary(Ω), joinpath(WPB.PROJECT_ROOT,"test","IO","sphere")) |> vtk_save
    rm(joinpath(@__DIR__,"sphere.vtu"))
end
