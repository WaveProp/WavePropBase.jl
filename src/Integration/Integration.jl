#=
Methods for integrating over instances of [`AbstractReferenceShape`](@ref).

Besides some standard quadrature rules used to integrate smooth functions, it also defines
singular integration routines useful for (weakly) singular integrands.
=#

include("quadrulestables.jl")
include("quadrule.jl")
include("adaptive.jl")
include("singularityhandler.jl")
include("singularquadrule.jl")
