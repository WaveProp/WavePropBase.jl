using Test
using WriteVTK
using Gmsh
import WavePropBase as WPB

@testset "VTK export" begin
    # This test should simply not throw an error
    WPB.clear_entities!()
    fname = joinpath(WPB.PROJECT_ROOT, "test", "IO", "sphere.msh")
    Ω, M = WPB.@gmsh begin
        WPB.gmsh_read_msh(fname)
    end
    vtk_save(WPB.vtk_mesh_file(M, joinpath(WPB.PROJECT_ROOT, "test", "IO", "ball")))
    rm(joinpath(WPB.PROJECT_ROOT, "test", "IO", "ball.vtu"))
    vtk_save(WPB.vtk_mesh_file(M, Ω, joinpath(WPB.PROJECT_ROOT, "test", "IO", "ball")))
    rm(joinpath(WPB.PROJECT_ROOT, "test", "IO", "ball.vtu"))
    vtk_save(WPB.vtk_mesh_file(M,
                               WPB.external_boundary(Ω),
                               joinpath(WPB.PROJECT_ROOT, "test", "IO", "sphere")))
    rm(joinpath(@__DIR__, "sphere.vtu"))
end

##
using StaticArrays
import WavePropBase as WPB
using WriteVTK
using Gmsh

WPB.clear_entities!()
gmsh.initialize()
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1)
gmsh.option.setNumber("Mesh.ElementOrder", 1)
gmsh.model.occ.synchronize()
gmsh.model.occ.addDisk(0,0,0,1,1)
gmsh.model.mesh.recombine()
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)
gmsh.model.mesh.recombine()
Ω = WPB.gmsh_import_domain(;dim=2)
M = WPB.gmsh_import_mesh(Ω;dim=2)
Γ = WPB.external_boundary(Ω)

WPB.vtk_mesh_file(M,Ω,"test") |> vtk_save

gmsh.finalize()

using CairoMakie
using Makie.GeometryBasics

els = M[collect(keys(M))[2]]
p = [Polygon(Point.(el.vals)) for el in els]
fig,axis,pl = poly(p)

els = M[collect(keys(M))[3]]
p = [Polygon(Point.(el.vals)) for el in els]
poly!(axis,p)

display(fig)


##

import WavePropBase as WPB
using StaticArrays
l1 = WPB.ParametricEntity(s-> SVector(s[1],0.1*sin(2π*s[1])),0,1)
l2 = WPB.ParametricEntity(s-> SVector(1+0.1*sin(2π*s[1]),s[1]),0,1)
l3 = WPB.ParametricEntity(s-> SVector(1-s[1], 1 + 0.1*sin(2π*s[1])),0,1)
l4 = WPB.ParametricEntity(s-> SVector(0.1*sin(2π*s[1]), 1-s[1]),0,1)
ent = WPB.transfinite_rectangle(l1,l2,l3,l4)

msh = WPB.meshgen(ent;meshsize=0.01)
els = msh[first(keys(msh))]

p = [Polygon(Point.(WPB.LagrangeElement(el).vals)) for el in els]
fig,axis,pl = poly(p)

el = els[1]



WPB.LagrangeElement(el)
pts = WPB.Point2D[]
lines = MeshCell[]
for k in keys(msh)
    for el in msh[k]
        push!(pts,el(0.0))
        push!(pts,el(1.0))
        cell = MeshCell(VTKCellTypes.VTK_LINE,[length(pts)-1,length(pts)])
        push!(lines,cell)
    end
end
npts = length(pts)
pts_mat = collect(reinterpret(reshape,Float64,pts))
vtk_grid("boundaries",pts_mat,lines) do vtk
    vtk["line_id"] = 1:length(lines)
end
