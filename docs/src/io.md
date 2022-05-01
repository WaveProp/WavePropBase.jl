# [IO module](@id io-section)

```@meta
CurrentModule = WavePropBase
```

## Overview

The `IO` module provides various *recipes* for `Plots.jl`, as well functionality
to export meshes and solutions to `.vtk` format using `WriteVTK`. Note that the
`vtkIO.jl` file is only loaded if `WriteVTK` is available (i.e. you must type
`using WriteVTK` for the file to be loaded).
