using WavePropBase
using WavePropBase.Geometry
using WavePropBase.IO
using WavePropBase.Mesh
using WriteVTK

@testset "VTK export" begin
    # This test should simply not throw an error
    Geometry.clear_entities!()
    Ω, M = WavePropBase.IO.gmsh_sphere()
    vtk_mesh_file(M, joinpath(@__DIR__,"ball")) |> vtk_save
    rm(joinpath(@__DIR__,"ball.vtu"))
    vtk_mesh_file(M, Ω, joinpath(@__DIR__,"ball")) |> vtk_save
    rm(joinpath(@__DIR__,"ball.vtu"))
    vtk_mesh_file(M, external_boundary(Ω), joinpath(@__DIR__,"sphere")) |> vtk_save
    rm(joinpath(@__DIR__,"sphere.vtu"))
end
