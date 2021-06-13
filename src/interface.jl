"""
    ambient_dimension(x)

Dimension of the ambient space where `x` lives. For geometrical objects this can
differ from its [`geometric_dimension`](@ref); for example a triangle in `ℝ³` has
ambient dimension `3` but geometric dimension `2`, while a curve in `ℝ³` has
ambient dimension 3 but geometric dimension 1.
"""
function ambient_dimension end
@interface ambient_dimension

"""
    geometric_dimension(x)
    geometric_dimension(Ω::Domain)

Number of degrees of freedom necessary to locally represent the geometrical
object. For example, lines have geometric dimension of 1 (whether in `ℝ²` or in
`ℝ³`), while surfaces have geometric dimension of 2.

When the argument is a `Domain`, return the largest geometric dimension
encoutered.
"""
function geometric_dimension end
@interface geometric_dimension

"""
    tag(::AbstractEntity)

Integer tag used to idetify geometrical entities.
"""
function tag end
@interface tag

"""
    boundary(ω)

Return the boundary of `ω`. For a mesh element gives the `d-1` dimensional
elements composing its boundary, while for an entity gives the corresponding
`d-1` dimensional entities.
"""
function boundary end
@interface boundary

"""
    diameter(Ω)

Largest distance between `x` and `y` for `x,y ∈ Ω`.
"""
function diameter end
@interface diameter

"""
    radius(Ω)

Half the [`diameter`](@ref).
"""
function radius end
@interface radius

"""
    center(Ω)
"""
function center end
@interface center

"""
    return_type(f)

The type returned by the function-like object `f`.
"""
function return_type end
@interface return_type

"""
    jacobian(F,x̂)

The Jacobian matrix `Aᵢⱼ = ∂Fᵢ/∂x̂ⱼ` at the parametric coordinate `x̂`.
"""
function jacobian end
@interface jacobian

"""
    normal(el,x̂)

The unit normal vector of `el` at the parametric coordinate `x̂`.
"""
function normal end
@interface normal

"""
    domain(f)

The domain of `f`. For elements of geometrical nature return the
`ReferenceShape` used to represent it.
"""
function domain end
@interface domain
