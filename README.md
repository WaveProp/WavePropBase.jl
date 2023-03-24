# WavePropBase

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://WaveProp.github.io/WavePropBase.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://WaveProp.github.io/WavePropBase.jl/dev)
[![Build
Status](https://github.com/WaveProp/WavePropBase/workflows/CI/badge.svg)](https://github.com/WaveProp/WavePropBase.jl/actions)
[![codecov](https://codecov.io/gh/WaveProp/WavePropBase.jl/branch/main/graph/badge.svg?token=codJo03vp6)](https://codecov.io/gh/WaveProp/WavePropBase.jl)

This package provides some basic functionality used across the
[`WaveProp`](https://github.com/WaveProp) organization. It defines a common set
of methods and structures to talk about things like domains, meshes, trees,
quadratures, etc. It is essentially a factorization of code that was being
duplicated across different packages. *Unless you are developing a package,
there is little reason for you to depend directly on `WavePropBase`.* See the
[documentation](https://WaveProp.github.io/WavePropBase.jl/stable) for more
details.
