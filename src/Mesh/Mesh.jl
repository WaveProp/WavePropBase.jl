#=

Some general mesh data structures used to solve PDEs and BIEs. Focus is on
generic representations and the possibility of handling curved meshes.

=#

include("abstractmesh.jl")
include("genericmesh.jl")
include("cartesianmesh.jl")
include("submesh.jl")
include("decompose.jl")
include("meshgen.jl")
include("gmsh.jl")
include("vtk.jl")
