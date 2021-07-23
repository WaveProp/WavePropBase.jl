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
    # methods
    reference_nodes,
    monomial_basis,
    lagrange_basis,
    barycentric_lagrange_matrix,
    barycentric_lagrange_weights,
    MonomialBasis,
    vandermond,
    lagrange_basis,
    PolynomialBasis,
    LagrangeBasis,
    gradLagrangeBasis,
    trace,
    grad

include("polynomials.jl")
include("tensorlaginterp.jl")
include("chebinterp.jl")
include("element.jl")
# include("monomial.jl")
# include("polynomialbasis.jl")

end # module Interpolation