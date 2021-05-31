# a place to document methods which are inteded to be extended. These are pushed
# into the global constant INTERFACE_LIST to facilitate processing such as e.g.
# importing all interface methods by other packages or generating a
# documentation page for this methods programatically

"""
    ambient_dimension(x)

Dimension of the ambient space where `x` lives. For geometrical objects this can
differ from its [`geometric_dimension`](@ref); for example a triangle in `ℝ³` has
ambient dimension `3` but geometric dimension `2`, while a curve in `ℝ³` has
ambient dimension 3 but geometric dimension 1.
"""
function ambient_dimension end
push!(INTERFACE_LIST,:ambient_dimension)

"""
    geometric_dimension(x)

Number of degrees of freedom necessary to locally represent the geometrical
object. For example, lines have geometric dimension of 1 (whether in `ℝ²` or in
`ℝ³`), while surfaces have geometric dimension of 2.
"""
function geometric_dimension end
push!(INTERFACE_LIST,:geometric_dimension)

"""
    boundary(ω)

Return the boundary of `ω`. For a mesh element gives the `d-1` dimensional
elements composing its boundary, while for an entity gives the corresponding
`d-1` dimensional entities.
"""
function boundary end
push!(INTERFACE_LIST,:boundary)
