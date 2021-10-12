"""
    module Integration

Methods for integrating over instances of [`AbstractReferenceShape`](@ref).

Besides some standard quadrature rules used to integrate smooth functions, it also defines
singular integration routines useful for (weakly) singular integrands.
"""
module Integration

using StaticArrays
using LinearAlgebra

using WavePropBase.Utils
using WavePropBase.Geometry
using WavePropBase.Interpolation

# import all methods in WavePropBase.INTERFACE_METHODS
using WavePropBase
WavePropBase.@import_interface

export
    # abstract types
    AbstractQuadratureRule,
    AbstractSingularityHandler,
    # types
    Gauss,
    Trapezoidal,
    TrapezoidalOpen,
    Fejer,
    TensorProductQuadrature,
    Kress,
    KressR,
    KressP,
    IMT,
    Duffy,
    SingularQuadratureRule,
    # functions
    integrate,
    qnodes,
    qweights,
    qnormals,
    integration_measure,
    qrule_for_reference_shape,
    singular_weights,
    singular_quadrature,
    TensorProductSingularityHandler

include("quadrulestables.jl")
include("quadrule.jl")
include("singularityhandler.jl")
include("singularquadrule.jl")

end
