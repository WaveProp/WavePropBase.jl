using Test
import WavePropBase as WPB
using CairoMakie
using Gmsh

WPB.clear_entities!()
Ω, msh = WPB.@gmsh begin
    WPB.clear_entities!()
    gmsh.model.occ.addDisk(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)
    Ω = WPB.gmsh_import_domain(; dim=2)
    M = WPB.gmsh_import_mesh(Ω; dim=2)
    Ω, M
end
Γ = WPB.boundary(Ω)

@test_nowarn lines(view(msh,Γ))
@test_nowarn poly(view(msh,Ω))

WPB.clear_entities!()
Ω, msh = WPB.@gmsh begin
    WPB.clear_entities!()
    gmsh.model.occ.addSphere(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    Ω = WPB.gmsh_import_domain(; dim=3)
    M = WPB.gmsh_import_mesh(Ω; dim=3)
    Ω, M
end
Γ = WPB.boundary(Ω)

@test_nowarn poly(view(msh,Γ),color=:lightgray,strokewidth=1)
@test_nowarn poly(view(msh,Ω),color=:lightgray,strokewidth=1,transparency=true)
