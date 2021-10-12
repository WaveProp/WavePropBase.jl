# [Interpolation module](@id interpolation-section)

```@meta
CurrentModule = WavePropBase
```

## Overview

The `Interpolation` module defines an interface for talking about polynomial
spaces and interpolation. The central concept of this module is that of an
[`AbstractElement{D,T}`](@ref Interpolation.AbstractElement), which maps points
on a reference domain (of type [`<:AbstractReferenceShape`](@ref
Geometry.AbstractReferenceShape)) to values of type `T`. The `AbstractElement`
interface expects the methods `(el::AbstractElement)(u)`, which evaluates the
underlying interpolant at the (reference) coordinate `u`, and `jacobian(el,u)`,
which computes the jacobian at the (reference) coordinate `u`.

`AbstractElement`s are commonly used to describe functions over
[`AbstractReferenceShape`](@ref Geometry.AbstractReferenceShape)s, and such
functions can in turn be used to describe more complex geometrical shapes used
e.g. in a mesh. By composing a function representation on a reference element
with the representation of the element itself as a map from the reference
element, we can therefore represent a function over a (possibly complex)
geometrical element.

The `Interpolation` module provides a few concrete implementations of
`AbstractElement`s which are described next, but other packages may use the
interface by implementing the methods `(::AbstractElement)(x)` and
`jacobian(el,x)` (see e.g. `ParametricElement` in the
[ParametricSurfaces](https://github.com/WaveProp/ParametricSurfaces.git)
package).

## LagrangeElements

One of the simplest `AbstractElement`s is the `LagrangeElement{D,T,Np}`, which
defines a polynomial mapping the `Np` reference nodes on `D` to `Np` values of
type `T`. The [`reference_nodes`](@ref Interpolation.reference_nodes) depend
only on `D` and `Np` (and therefore on the type of the element). We use the same
convention as *Gmsh* to define the order of the reference nodes on the various
reference shapes; see [node
ordering](https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering) on the *Gmsh*
documentation for a more in-depth description. 

In `WaveProp`, `LagrangeElement`s are often used to describe (possibly curved)
mesh elements. The triangle with vertices `(1,1),(2,2),(1.5,3)` can for example
be created using:

```@example triangle-element
using WavePropBase, Plots
pyplot() # hide
pts = (1,1),(2,2),(1.5,3)
el  = Interpolation.LagrangeTriangle(pts)
plot(el)
```

Note that, as per the `AbstractElement` interface, you may evaluate the
parametrization and the jacobian at any point on the reference element (not just
the reference nodes):

```@example triangle-element
u = Geometry.Point2D(0.25,0.25)
@show el(u), jacobian(el,u)
```

If the triangle is instead in three dimensions, it suffices to pass
three-dimensional points:

```@example triangle-element
pts = (1,1,0),(2,2,1),(1.5,3,1)
el  = Interpolation.LagrangeTriangle(pts)
plot(el)
```

Very similar constructs can be used to work higher order (curved) triangles, or
with other `LagrangeElements` such as [`LagrangeSquare`](@ref
Interpolation.LagrangeSquare) or [`LagrangeTetrahedron`](@ref
Interpolation.LagrangeTetrahedron); see their docstrings for more details.

```@index
Modules = [WavePropBase.Interpolation]
```
