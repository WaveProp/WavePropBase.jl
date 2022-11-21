#=
Tabulated Gaussian quadrature rules from GMSH.
Obtained from `https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_7_0/Numeric/GaussQuadratureTri.cpp`
and `https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_7_0/Numeric/GaussQuadratureTet.cpp`.
=#

include("quadrulestables_tetrahedron.jl")
include("quadrulestables_triangle.jl")

# tabulated G7K15 rules
const SEGMENT_G7 = (SVector{7}(SVector(0.025446043828620757),
                               SVector(0.12923440720030277),
                               SVector(0.2970774243113014),
                               SVector(0.5),
                               SVector(0.7029225756886985),
                               SVector(0.8707655927996972),
                               SVector(0.9745539561713792)),
                    SVector{7}(0.06474248308443491,
                               0.13985269574463835,
                               0.19091502525255943,
                               0.20897959183673454,
                               0.19091502525255943,
                               0.13985269574463835,
                               0.06474248308443491))

# Kronrod 15 point rule
const SEGMENT_K15 = (SVector{15}(SVector(0.025446043828620757),
                                 SVector(0.12923440720030277),
                                 SVector(0.2970774243113014),
                                 SVector(0.5),
                                 SVector(0.7029225756886985),
                                 SVector(0.8707655927996972),
                                 SVector(0.9745539561713792),
                                 SVector(0.004272314439593694),
                                 SVector(0.06756778832011551),
                                 SVector(0.20695638226615443),
                                 SVector(0.39610752249605075),
                                 SVector(0.6038924775039493),
                                 SVector(0.7930436177338456),
                                 SVector(0.9324322116798844),
                                 SVector(0.9957276855604063)),
                     SVector{15}(0.03154604631498921,
                                 0.07032662985776296,
                                 0.09517528903239279,
                                 0.10474107054236396,
                                 0.09517528903239279,
                                 0.07032662985776296,
                                 0.03154604631498921,
                                 0.011467661005264628,
                                 0.052395005161125115,
                                 0.08450236331963394,
                                 0.10221647003764939,
                                 0.10221647003764939,
                                 0.08450236331963394,
                                 0.052395005161125115,
                                 0.011467661005264628))

# Radon' 7 point rule of order 5 for triangle
const TRIANGLE_R5N7 = (SVector{7}(SVector(0.33333333333333333, 0.33333333333333333),
                                  SVector(0.79742698535308720, 0.10128650732345633),
                                  SVector(0.10128650732345633, 0.79742698535308720),
                                  SVector(0.10128650732345633, 0.10128650732345633),
                                  SVector(0.05971587178976981, 0.47014206410511505),
                                  SVector(0.47014206410511505, 0.05971587178976981),
                                  SVector(0.47014206410511505, 0.47014206410511505)),
                       SVector{7}(0.22500000000000000 / 2,
                                  0.12593918054482717 / 2,
                                  0.12593918054482717 / 2,
                                  0.12593918054482717 / 2,
                                  0.13239415278850616 / 2,
                                  0.13239415278850616 / 2,
                                  0.13239415278850616 / 2))

# Laurie's 19 point rule of order 8
const TRIANGLE_L8N19 = (SVector{19}(SVector(0.3333333333333333, 0.3333333333333333),
                                    SVector(0.7974269853530872, 0.1012865073234563),
                                    SVector(0.1012865073234563, 0.7974269853530872),
                                    SVector(0.1012865073234563, 0.1012865073234563),
                                    SVector(0.0597158717897698, 0.4701420641051151),
                                    SVector(0.4701420641051151, 0.0597158717897698),
                                    SVector(0.4701420641051151, 0.4701420641051151),
                                    SVector(0.5357953464498992, 0.2321023267750504),
                                    SVector(0.2321023267750504, 0.5357953464498992),
                                    SVector(0.2321023267750504, 0.2321023267750504),
                                    SVector(0.9410382782311209, 0.0294808608844396),
                                    SVector(0.0294808608844396, 0.9410382782311209),
                                    SVector(0.0294808608844396, 0.0294808608844396),
                                    SVector(0.7384168123405100, 0.2321023267750504),
                                    SVector(0.7384168123405100, 0.0294808608844396),
                                    SVector(0.2321023267750504, 0.7384168123405100),
                                    SVector(0.2321023267750504, 0.0294808608844396),
                                    SVector(0.0294808608844396, 0.7384168123405100),
                                    SVector(0.0294808608844396, 0.2321023267750504)),
                        SVector{19}(0.0378610912003147,
                                    0.0376204254131829,
                                    0.0376204254131829,
                                    0.0376204254131829,
                                    0.0783573522441174,
                                    0.0783573522441174,
                                    0.0783573522441174,
                                    0.1162714796569659,
                                    0.1162714796569659,
                                    0.1162714796569659,
                                    0.0134442673751655,
                                    0.0134442673751655,
                                    0.0134442673751655,
                                    0.0375097224552317,
                                    0.0375097224552317,
                                    0.0375097224552317,
                                    0.0375097224552317,
                                    0.0375097224552317,
                                    0.0375097224552317) / 2)

##
# Dictionaries that contains the quadrature rules
# for various `AbstractReferenceShape`. The dictionary
# key represents the number of quadrature nodes.

const TRIANGLE_GAUSS_QRULES = Dict(1 => TRIANGLE_G1N1,
                                   3 => TRIANGLE_G2N3,
                                   4 => TRIANGLE_G3N4,
                                   6 => TRIANGLE_G4N6,
                                   7 => TRIANGLE_G5N7,
                                   12 => TRIANGLE_G6N12,
                                   13 => TRIANGLE_G7N13,
                                   16 => TRIANGLE_G8N16,
                                   19 => TRIANGLE_G9N19,
                                   25 => TRIANGLE_G10N25,
                                   33 => TRIANGLE_G12N33)

# map a desired quadrature order to the number of nodes
const TRIANGLE_GAUSS_ORDER_TO_NPTS = Dict(1 => 1,
                                          2 => 3,
                                          3 => 4,
                                          4 => 6,
                                          5 => 7,
                                          6 => 12,
                                          7 => 13,
                                          8 => 16,
                                          9 => 19,
                                          10 => 25,
                                          12 => 33)
const TRIANGLE_GAUSS_NPTS_TO_ORDER = Dict((v, k) for (k, v) in TRIANGLE_GAUSS_ORDER_TO_NPTS)

const TETAHEDRON_GAUSS_QRULES = Dict(1 => TETAHEDRON_G1N1,
                                     4 => TETAHEDRON_G2N4,
                                     5 => TETAHEDRON_G3N5,
                                     11 => TETAHEDRON_G4N11,
                                     14 => TETAHEDRON_G5N14,
                                     24 => TETAHEDRON_G6N24,
                                     31 => TETAHEDRON_G7N31,
                                     43 => TETAHEDRON_G8N43)

# map a desired quadrature order to the number of nodes
const TETRAHEDRON_GAUSS_ORDER_TO_NPTS = Dict(1 => 1, 2 => 4, 3 => 5, 4 => 11, 5 => 14,
                                             6 => 24, 7 => 31, 8 => 43)
const TETRAHEDRON_GAUSS_NPTS_TO_ORDER = Dict((v, k)
                                             for (k, v) in TETRAHEDRON_GAUSS_ORDER_TO_NPTS)

##
const GAUSS_QRULES = Dict(ReferenceTriangle => TRIANGLE_GAUSS_QRULES,
                          ReferenceTetrahedron => TETAHEDRON_GAUSS_QRULES)

"""
    _get_gauss_and_qweights(R::Type{<:AbstractReferenceShape{D}}, N) where D

Returns the `N`-point symmetric gaussian qnodes and qweights `(x, w)` for integration over `R`.
"""
function _get_gauss_qcoords_and_qweights(R::Type{<:AbstractReferenceShape{D}}, N) where {D}
    if !haskey(GAUSS_QRULES, R) || !haskey(GAUSS_QRULES[R], N)
        error("quadrature rule not found")
    end
    qrule = GAUSS_QRULES[R][N]
    @assert length(qrule) == N
    # qnodes
    qnodestype = SVector{N,SVector{D,Float64}}
    x = qnodestype([q[1] for q in qrule])
    # qweights
    qweightstype = SVector{N,Float64}
    w = qweightstype([q[2] for q in qrule])
    return x, w
end
