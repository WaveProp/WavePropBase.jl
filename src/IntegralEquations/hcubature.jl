#=
Routines for correcting the singular or nearly singular entries of an integral
operator based on:
    1. Identifying apriorily the entries of the integral operator which need to
    be corrected. This is done based on distance between points/elements and
    analytical knowldge of the convergence rate of the quadrature and the
    underlying kernels
    2. Performing a change-of-variables to alleviate the singularity of the
    integrand
    3. Doing an adaptive quadrature on the new integrand using HCubature.jl

The specialized integration is precomputed on the reference element for each
quadrature node of the standard quadrature.
=#

function hcubature_correction(iop::IntegralOperator, shandler_dict=nothing; max_dist,
                              hcubature_kwargs...)
    X, Y = target_mesh(iop), source_mesh(iop)
    Xqnodes = qnodes(X)
    Yqnodes = qnodes(Y)
    K = kernel(iop)
    T = eltype(iop)
    dict_near = near_interaction_list(Xqnodes, Y; tol=max_dist)
    Is, Js, Vs = Int[], Int[], T[]
    for E in keys(Y)
        N = geometric_dimension(E)
        a, b = ntuple(i -> 0, N), ntuple(i -> 1, N)
        τ̂ = domain(E)
        iter = Y[E] # iterator over elements of type E in the source mesh
        qreg = etype2qrule(Y, E)
        shand = isnothing(shandler_dict) ? identity : shandler_dict[E]
        L = lagrange_basis(qreg)
        x̂, ŵ = qreg()
        el2qtags = etype2qtags(Y, E)
        nearlist = dict_near[E]
        buffer = hcubature_buffer(x -> one(T) * L(a) * first(ŵ), a, b)
        # TODO: add function barrier to make the loop below type-stable
        for (n, el) in enumerate(iter)
            jglob = view(el2qtags, :, n)
            for i in nearlist[n]
                xnode = Xqnodes[i]
                # closest quadrature node
                dmin, j = findmin(n -> norm(coords(xnode) - coords(Yqnodes[jglob[n]])),
                                  1:length(jglob))
                x̂nearest = x̂[j]
                dmin > max_dist && continue
                # FIXME: better estimate the distance above using a bounding sphere on
                # `el` instead of the smallest distance to the quadrature nodes
                W = sum(decompose(τ̂, x̂nearest)) do l
                    I, E = hcubature(a, b; buffer, hcubature_kwargs...) do ŷs
                        ŷ = l(ŷs)
                        l′ = integration_measure(l, ŷs)
                        y = el(ŷ)
                        jac = jacobian(el, ŷ)
                        ν = _normal(jac)
                        τ′ = _integration_measure(jac)
                        return K(xnode, (coords=y, normal=ν)) * L(ŷ) * τ′ * l′
                    end
                    return I
                end
                for (k, j) in enumerate(jglob)
                    push!(Is, i)
                    push!(Js, j)
                    push!(Vs, W[k] - iop[i, j])
                end
            end
        end
    end
    return sparse(Is, Js, Vs)
end
