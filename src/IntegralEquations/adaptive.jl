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

function singularquadrule_correction(iop::IntegralOperator, qsing_dict; tol)
    X, Y = target_mesh(iop), source_mesh(iop)
    Xqnodes = qnodes(X)
    Yqnodes = qnodes(Y)
    K = kernel(iop)
    T = eltype(iop)
    dict_near = near_interaction_list(Xqnodes, Y; tol)
    Is, Js, Vs = Int[], Int[], T[]
    for E in keys(Y)
        τ̂ = domain(E)
        iter = Y[E]
        qreg = etype2qrule(Y, E)
        qsing = qsing_dict[E]
        x̂,ŵ   = qreg()
        el2qtags = etype2qtags(Y, E)
        nearlist = dict_near[E]
        L = lagrange_basis(qreg)
        for (n, el) in enumerate(iter)
            jglob = view(el2qtags, :, n)
            for i in nearlist[n]
                xnode = Xqnodes[i]
                # closest quadrature node
                dmin,j = findmin(n -> norm(coords(xnode) - coords(Yqnodes[jglob[n]])),
                           1:length(jglob))
                x̂nearest = x̂[j]
                dmin > tol && continue # FIXME: better estimate the distance using a bounding sphere
                W = sum(decompose(τ̂,x̂nearest)) do l
                    integrate(qsing) do ŷs
                        ŷ   = l(ŷs)
                        l′  = integration_measure(l,ŷs)
                        y   = el(ŷ)
                        jac = jacobian(el, ŷ)
                        ν = _normal(jac)
                        τ′ = _integration_measure(jac)
                        any(isnan,L(ŷ)) && (error("Nan found $(ŷ), $(L(ŷ))"))
                        K(xnode,(coords=y,normal=ν))*L(ŷ)*τ′*l′
                    end
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
