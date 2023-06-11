# [Geometry](@id geometry-section)

```@meta
CurrentModule = WavePropBase
```

## Overview

Several backends for the actual representation of geometrical objects are
supported by `WavePropBase`. The most basic backend is the
[`ParametricEntity`](@ref) type, describes geometrica entities for which an
explicit (user-defined) parametrization is available. This is useful for
describing simple shapes such as circles, ellipses, and so on. Dealing with more complex
shapes require a proper CAD engine, and `WavePropBase` does not possess such
features. You can, however, use the `Gmsh` to create an import both geometrical
entities and meshes---see the [Gmsh](@ref gmsh-section) for more information.

```@example
using Plots
import WavePropBase as WPB
ent = WPB.ParametricEntity(s->(cos(s[1]),sin(s[1])),0,2π) 
ss  = range(0,2π,length=100)
pts = ent.(ss)
fig = plot([pt[1] for pt in pts], [pt[2] for pt in
pts],aspect_ratio=:equal,label="")
```

```@autodocs
Modules = [WavePropBase]
Pages   = ["Geometry/domain.jl","Geometry/entities.jl","Geometry/point.jl","Geometry/referenceshapes.jl"]
```