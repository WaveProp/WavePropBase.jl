#=
Singular quadrarute technique based on change of variables to compute a sparse correction to integral
operators when the kernel is weakly singular. The technique consists of:
    1. Identifying a priorily the entries of the integral operator which need
    correction. Those are given by entries whose target point and source element
    lie closer than a certain threshold
    2. For each element in the source surface, and for each target point close
    to that element, compute the matrix entries by performing an specialized
    integration on the Lagrange interpolants.

The specialized integration is precomputed on the reference element for each
quadrature node of the standard quadrature.
=#

function singularquadrule_correction(iop::IntegralOperator, qsing; tol)
    X, Y = target_mesh(iop), source_mesh(iop)
    Xqnodes = qcoords(X)
    Yqnodes = qcoords(Y)
    K = kernel(iop)
    T = eltype(iop)
    dict_near = near_interaction_list(Xqnodes, Y; tol)
    Is, Js, Vs = Int[], Int[], T[]
    for E in keys(Y)
        iter = Y[E]
        qreg = etype2qrule(Y, E)
        # qsing = qsing_dict[E]
        el2qtags = etype2qtags(Y, E)
        nearlist = dict_near[E]
        nodes, weights = singular_quadrature_suite(qreg, qsing)
        L = lagrange_basis(qreg)
        for (n, el) in enumerate(iter)
            jglob = view(el2qtags, :, n)
            for i in nearlist[n]
                xnode = Xqnodes[i]
                dmin,j = findmin(n -> norm(coords(xnode) - coords(Yqnodes[jglob[n]])),
                           1:length(jglob))
                dmin > tol && continue
                W = sum(zip(nodes[j], weights[j])) do (ŷ, ŵ)
                    y = el(ŷ)
                    jac = jacobian(el, ŷ)
                    ν = _normal(jac)
                    μ = _integration_measure(jac)
                    K(xnode,(coords=y,normal=ν))*ŵ*μ*L(ŷ)
                end
                for (k,j) in enumerate(jglob)
                    push!(Is, i)
                    push!(Js, j)
                    push!(Vs, W[k] - iop[i,j])
                end
            end
        end
    end
    return sparse(Is, Js, Vs)
end

# #=
# Adaptive quadrature technique to compute a sparse correction to integral
# operators when the kernel is weakly singular. The technique consists of:
#     1. Identifying a priorily the entries of the integral operator which need
#     correction. Those are given by entries whose target point and source element
#     lie closer than a certain threshold
#     2. For each element in the source surface, and for each target point close
#     to that element, compute the matrix entries by performing an adaptive
#     integration
# =#

# """
#     near_interaction_list_adaptive(X,Y::AbstractMesh; atol)

# For each element `el` of type `E` in `Y`, return the indices of the points in `X` which
# are closer than `atol` to the `center` of `el`.

# This function returns a dictionary where e.g. `Dict[E][5] --> Vector{Int}` gives
# the indices of points in `X` which are closer than atol to the center of the
# fifth element of type `E`.
# """
# function near_interaction_list(X, M::AbstractMesh; atol)
#     dict = Dict{DataType,Vector{Vector{Int}}}()
#     for E in keys(M)
#         pts = [center(el) for el in M[E]]
#         kdtree = KDTree(pts)
#         push!(dict, E => inrange(kdtree, X, atol))
#     end
#     return dict
# end

# function adaptive_correction(iop::IntegralOperator)
#     X, Y = target_surface(iop), source_surface(iop)
#     # maximum element size
#     lmax = -Inf
#     dict_near = near_interaction_list(dofs(X), Y; atol=0.1)
#     return
# end

# function assemble_gk(iop; compress=Matrix)
#     X, Y = target_surface(iop), source_surface(iop)
#     @timeit_debug "assemble dense part" begin
#         out = compress(iop)
#     end
#     @timeit_debug "compute near interaction list" begin
#         dict_near = near_interaction_list(dofs(X), Y; atol=0.1)
#     end
#     @timeit_debug "compute singular correction" begin
#         correction = singular_weights_gk(iop, dict_near)
#     end
#     return axpy!(1, correction, out) # out <-- out + correction
# end

# function singular_weights_gk(iop::IntegralOperator, dict_near)
#     X, Y = target_surface(iop), source_surface(iop)
#     T = eltype(iop)
#     Is = Int[]
#     Js = Int[]
#     Vs = T[]
#     # loop over all integration elements by types
#     for E in keys(Y)
#         iter = ElementIterator(Y, E)
#         qreg = Y.etype2qrule[E]
#         lag_basis = lagrange_basis(qnodes(qreg))
#         @timeit_debug "singular weights" begin
#             _singular_weights_gk!(Is, Js, Vs, iop, iter, dict_near, lag_basis, qreg)
#         end
#     end
#     Sp = sparse(Is, Js, Vs, size(iop)...)
#     return Sp
# end

# function _singular_weights_gk!(Is, Js, Vs, iop, iter, dict_near, lag_basis, qreg)
#     yi = qnodes(qreg)
#     K = kernel(iop)
#     X = target_surface(iop)
#     Y = source_surface(iop)
#     E = eltype(iter)
#     list_near = dict_near[E]
#     # Over the lines below we loop over all element, find which target nodes
#     # need to be regularized, and compute the regularized weights on the source
#     # points for that target point.
#     for (n, el) in enumerate(iter)               # loop over elements
#         jglob = Y.elt2dof[E][:, n]
#         for (i, jloc) in list_near[n]            # loop over near targets
#             ys = yi[jloc]
#             xdof = dofs(X)[i]
#             for (m, l) in enumerate(lag_basis)   # loop over lagrange basis
#                 val, _ = quadgk(0, ys[1], 1; atol=1e-14) do ŷ
#                     y = el(ŷ)
#                     jac = jacobian(el, ŷ)
#                     μ = integration_measure(jac)
#                     ydof = NystromDOF(y, -1.0, jac, -1, -1)
#                     return K(xdof, ydof) * μ * l(ŷ)
#                 end
#                 # push data to sparse matrix
#                 push!(Is, i)
#                 j = jglob[m]
#                 push!(Js, j)
#                 push!(Vs, val - iop[i, j])
#             end
#         end
#     end
#     return Is, Js, Vs
# end
