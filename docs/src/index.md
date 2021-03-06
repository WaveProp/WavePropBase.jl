
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

## Package structure

The package consists of a single module, and has no exported structures. The
code has been physically organized into the following subfolders:

- [Utils](@ref utils-section)
- [Geometry](@ref geometry-section)
- [Trees](@ref trees-section)
- [Interpolation](@ref interpolation-section)
- [Integration](@ref integration-section)
- [Mesh](@ref mesh-section)
- [IO](@ref io-section)

## Index

```@index
Modules = [WavePropBase]
```
