# [Integration module](@id integration-section)

## Overview

The `Integration` module provides various quadrature rules and routines for
integrating functions defined over [`AbstractReferenceShape`](@ref
Geometry.AbstractReferenceShape)s. Quadrature rules are expected to inherit from
the `AbstractQuadratureRule{D}`, which describes a set of nodes and weights used
to integrate a function over a reference domain `D<:AbstractReferenceShape`.
Subtypes of `AbstractQuadratureRule` should implement at least
`(::AbstractQuadratureRule)()` returning the nodes and weights.

The
`Integration` module defines various quadrature rules such as

```@index
Modules = [WavePropBase.Integration]
```

## Regular integration rules

TODO

## Singular integration

TODO
