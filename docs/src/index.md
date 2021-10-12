
# WavePropBase

```@meta
CurrentModule = WavePropBase
```

## Overview

This package provides some basic functionality used across the
[`WaveProp`](https://github.com/WaveProp) organization. It defines a common set
of methods and structures to talk about things like domains, meshes, trees,
quadratures, etc. It is essentially a factorization of code that was being
duplicated across different packages. *Unless you are developing a package, there
is little reason for you to depend directly on `WavePropBase`.*

## Interface

The methods which are to be extended by other packages are defined at
the top-level module (`WavePropBase`), and are grouped for convenience in the
global variable [`INTERFACE`](@ref). You can automatically import/export all
the interface methods using [`@import_interface`](@ref) and
[`@export_interface`](@ref).

## Submodules

Most of the functionality has been organized inside the following sub-modules:

- [Utils](@ref utils-section)
- [Geometry](@ref geometry-section)
- [Trees](@ref trees-section)
- [Interpolation](@ref interpolation-section)
- [Integration](@ref integration-section)
- [Mesh](@ref mesh-section)
- [IO](@ref io-section)

Refer to each module's documentation for further information on how to use them,
as well as to an index of their methods/structures.

## Index

```@index
Modules = [WavePropBase]
```
