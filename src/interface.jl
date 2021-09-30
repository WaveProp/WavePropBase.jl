const INTERFACE = Vector{Symbol}()

# Maybe just write the text into a file instead of using a macro here?
macro import_interface()
    ex = Expr(:block)
    ex.args = [:(import WavePropBase: $f) for f in INTERFACE]
    return ex
end

macro export_interface()
    ex = Expr(:block)
    ex.args = [:(export $f) for f in INTERFACE]
    return ex
end

"""
    ambient_dimension(x)

Dimension of the ambient space where `x` lives. For geometrical objects this can
differ from its [`geometric_dimension`](@ref); for example a triangle in `ℝ³` has
ambient dimension `3` but geometric dimension `2`, while a curve in `ℝ³` has
ambient dimension 3 but geometric dimension 1.
"""
function ambient_dimension end
push!(INTERFACE,:ambient_dimension)

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
push!(INTERFACE,:geometric_dimension)

"""
    dimension(space)

The length of a basis for `space`; i.e. the number of linearly independent elements
required to span `space`.
"""
function dimension end
push!(INTERFACE,:dimension)

"""
    boundary(ω)

Return the boundary of `ω`. For a mesh element gives the `d-1` dimensional
elements composing its boundary, while for an entity gives the corresponding
`d-1` dimensional entities.
"""
function boundary end
push!(INTERFACE,:boundary)

"""
    diameter(Ω)

Largest distance between `x` and `y` for `x,y ∈ Ω`.
"""
function diameter end
push!(INTERFACE,:diameter)

"""
    distance(Ω₁,Ω₂)

Minimal Euclidean distance between `Ω₁` and `Ω₂`.
"""
function distance end
push!(INTERFACE,:distance)

"""
    radius(Ω)

Half the [`diameter`](@ref).
"""
function radius end
push!(INTERFACE,:radius)

"""
    center(Ω)

Center of the smallest possible ball containing `Ω`.
"""
function center end
push!(INTERFACE,:center)

"""
    return_type(f)

The type returned by function-like objects.
"""
function return_type end
push!(INTERFACE,:return_type)

"""
    jacobian(f,x)

Given a (possibly vector-valued) function `f : 𝐑ᵐ → 𝐅ᵐ`, return the `m × n`
matrix `Aᵢⱼ = ∂fᵢ/∂xⱼ`.
"""
function jacobian end
push!(INTERFACE,:jacobian)

"""
    normal(el,x̂)
    normal(jac::SMatrix)

The unit normal vector at coordinate `x̂`, guaranteed to be orthogonal to all
columns of `jacobian(el,x)`.
"""
function normal end
push!(INTERFACE,:normal)

"""
    domain(f)

The domain of the function `f`; i.e. if `f: Ω → R`, return `Ω`.
"""
function domain end
push!(INTERFACE,:domain)

"""
    image(f)

The image of the function `f`; i.e. if `f: Ω → R`, return `R`.
"""
function image end
push!(INTERFACE,:image)

"""
    parametrization(el)

Return the underlying parametrization of `el`.
"""
function parametrization end
push!(INTERFACE,:parametrization)

"""
    coords(x)

Return an `SVector` giving a cartesian coordinate used to sort `x`. For points,
this is simply the coordinates of the point, while for elements this may be the
center of the element. You must overload this function for your own type if you
want e.g. the clustering algorithms to work for `x`.
"""
function coords end
push!(INTERFACE,:coords)
