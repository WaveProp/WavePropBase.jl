# TODO: open an issue to discuss how this module should be implemented and what
# functionalities it should provide.
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
    cheb1weights,
    grad,
    vals

include("polynomials.jl")
include("tensorlaginterp.jl")
include("element.jl")

end # module Interpolation
