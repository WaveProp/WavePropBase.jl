# [Geometry](@id geometry-section)

```@meta
CurrentModule = WavePropBase
```
## Overview

The `Geometry` module is of fundamental importance as it defines the interface
expected from geometrical objects. Appropriately defining the geometry is one of
the first stages in setting up a simulation, and it precedes the creation of a
mesh. Recall that (at a high-level), the usual workflow is:

```math
    \fbox{Geometry} \rightarrow \fbox{Mesh} \rightarrow \fbox{Solver} \rightarrow \fbox{Solution}
```

`Geometry` then handles everything related to the representation of geometrical
entities and domains, as well as various operations on them. Furthermore, the
`Geometry` module also defines various simple reference shapes (see.
[`AbstractReferenceShape`](@ref)) which form the basis for interpolation and
integration procedures on more complex surfaces.

The most elementary object is an **entity**, which describes an atomic piece of
curve/surface/volume. Entities can then be grouped together to form **domains**,
and certain set operations can be performed on domains. Domains can later be
discretized to generate a mesh where actual computations can be performed. In
the following sections we will briefly describe the main structures of the
`Geometry` module, as well as provide some examples of creating geometrical
objects.

!!! note
    The `Geometry` module itself does not handle the actual representation of
    geometrical objects. Such functionality, which can be quite complex in
    practice, is delegated to other packages and/or software which implement the
    expected interface. Currently,
    [`ParametricSurfaces`](https://github.com/WaveProp/ParametricSurfaces) and
    [`GmshSDK`](https://github.com/WaveProp/GmshSDK) provide two ways of
    actually creating geometrical objects, and some examples in this page will
    use them to illustrate the ideas.

## Entities

Entities are the most basic geometrical objects. They can represent a point,
curve, surface, or volume. All entities are expected to inherit from
[`AbstractEntity`](@ref), and should extend the [`key`](@ref) and the
[`boundary`](@ref) methods. Calling `key(ent::AbstractEntity)` should return a
`(dim::Integer,tag::Integer)` which uniquely identifies the entity. The `dim`
value is simply its geometrical dimension: zero for points, one for curves, two
for surfaces, and so on, while `tag` is an integer used to distinguish the
entity from others.

Because the `(dim,tag)` key should be unique, the global variable `TAGS` exists
to keep track of the existing tags, and all types which inherit from
`AbstractEntity` should implement an inner constructor which calls
[`global_add_entity!`](@ref) upon creation of a new object. This function will
add the entity to the global dictionary [`ENTITIES`](@ref) so that it can always
be retrieved by its `(dim,tag)` key, as well as update the `TAGS` variable with
the new `(dim,tag)`. When running code in the `REPL`, it is sometimes useful to
call [`clear_entities!()`](@ref) to empty the `TAGS` and `ENTITIES` variables in
order to avoid unnecessary cluttering.

The `Geometry` module provides one (minimal) implementation of `AbstractEntity`:
the [`ElementaryEntity`](@ref) type. This type is useful for geometrical objects
for which no parametric information is available. For instance, when reading
files from *Gmsh*, the *Gmsh* entities are imported as `ElementaryEntity`s
containing a `dim`, `tag`, and `boundary` field, but no information on the
entity's parametrization is explicitly available. 

!!! tip
    The *entity* concept is widely used in *Gmsh* to describe the elementary
    geometrical objects inside a model, and the `ElementaryEntity` structure
    closely mimics what you may recover from the [Gmsh
    API](https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-API).    

## Domains

The [`Domain`](@ref) type provides a convenient way to groups several
`AbstractEntity` in order to form more complex geometrical objects. `Domain`s
support various set operations such as unions, intersection, and set difference.
Furthremore, domains can be used to index parts of a mesh, as described in the
[Mesh section](@ref mesh-section).

## Reference shapes

The `AbstractReferenceShape{N}` type describes singleton types representing
(fixed) geometrical shapes in `N` dimensions. Concrete subtypes include
[`ReferenceLine`](@ref), [`ReferenceTriangle`](@ref), [`ReferenceSquare`](@ref),
and the [`ReferenceTetrahedron`](@ref).

Various interpolation and integration routines can then be efficiently defined
on these singleton types, and more complex interpolation/integration over more
complex elements can be carried out by combining the routines on the reference
element with a parametrization of the complex element.
