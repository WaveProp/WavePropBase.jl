"""
    module Interpolation

Module implementing various interpolation schemes commonly used in finite
element and boundary element methods.
"""
module Interpolation

using StaticArrays
using LinearAlgebra

using WavePropBase
using WavePropBase.Utils
using WavePropBase.Geometry

WavePropBase.@import_interface

export
    # abstract types
    AbstractElement,
    AbstractHyperRectangle,
    # structs
    ParametricElement,
    LagrangeElement,
    LagrangePoint,
    LagrangeLine,
    LagrangeTriangle,
    LagrangeSquare,
    LagrangeTetrahedron,
    Pk,
    TensorLagInterp,
    ChebInterp,
    MonomialBasis,
    PolynomialBasis,
    LagrangeBasis,
    HyperRectangle,
    HyperCube,
    # methods
    reference_nodes,
    monomial_basis,
    lagrange_basis,
    barycentric_lagrange_matrix,
    barycentric_lagrange_weights,
    vandermond,
    lagrange_basis,
    gradLagrangeBasis,
    interpolation_nodes,
    trace,
    cheb1nodes,
    cheb1nodes_iter,
    cheb1weights,
    cheb2nodes,
    cheb2nodes_iter,
    cheb2weights,
    grad,
    vals,
    low_corner,
    high_corner,
    precompute_nodes_and_weights!

include("polynomials.jl")
include("element.jl")
include("hyperrectangle.jl")
include("tensorlaginterp.jl")

end # module Interpolation
