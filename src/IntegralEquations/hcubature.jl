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

"""
    hcubature_correction(iop::IntegralOperator; max_dist, kwargs...)

Given an integral operator `iop`, this function provides a sparse correction to
`iop` for the entries `i,j` such that the distance between the `i`-th target and
the `j`-th source is less than `max_dist`.

The correction is computed by using an adaptive quadrature based on the
`hcubature` routine from `HCubature.jl`; the `kwargs` are passed to `hcubature`,
and (see [its documentation](https://github.com/JuliaMath/HCubature.jl) for more
details).
"""
function hcubature_correction(iop::IntegralOperator; max_dist, kwargs...)
    # unpack type-unstable fields in iop, allocate output, and dispatch
    X, Y, K = target_mesh(iop), source_mesh(iop), kernel(iop)
    dict_near = near_interaction_list(X, Y; tol=max_dist)
    T  = eltype(iop)
    out = (I=Int[], J=Int[], V=T[])
    for E in keys(Y)
        # dispatch on element type
        nearlist = dict_near[E]
        qreg = etype2qrule(Y, E)
        L = lagrange_basis(qreg)
        iter = Y[E]
        _hcubature_correction_etype!(out, iter, qreg, L, nearlist, X, Y, K, max_dist, kwargs)
    end
    return sparse(out...)
end

@noinline function _hcubature_correction_etype!(out,iter,qreg,L,nearlist,X,Y,K,max_dist,kwargs)
    atol = haskey(kwargs,:atol) ? kwargs[:atol] : Inf
    E = eltype(iter)
    Xqnodes = qnodes(X)
    Yqnodes = qnodes(Y)
    N = geometric_dimension(E)
    a, b = ntuple(i -> 0, N), ntuple(i -> 1, N)
    τ̂ = domain(E)
    x̂, ŵ = collect.(qreg())
    el2qtags = etype2qtags(Y, E)
    buffer = hcubature_buffer(x -> one(eltype(out.V)) * L(a) * first(ŵ), a, b)
    for (n, el) in enumerate(iter)
        jglob = view(el2qtags, :, n)
        inear = nearlist[n]
        for i in inear
            xnode = Xqnodes[i]
            # closest quadrature node
            dmin, j = findmin(n -> norm(coords(xnode) - coords(Yqnodes[jglob[n]])),1:length(jglob))
            x̂nearest = x̂[j]
            dmin > max_dist && continue
            # FIXME: better estimate the distance above using a bounding sphere on
            # `el` instead of the smallest distance to the quadrature nodes
            μ = N == 1 ? Kress{2}() : N == 2 ? Duffy() : nothing
            # μ = identity
            ll = decompose(τ̂, x̂nearest)
            W = sum(ll) do l
                val, er = hcubature(a, b; buffer, kwargs...) do ŷs
                    ŷ = l(μ(ŷs))
                    l′ = integration_measure(l, μ(ŷs))
                    μ′ = integration_measure(μ,ŷs)
                    y  = el(ŷ)
                    jac = jacobian(el, ŷ)
                    ν = _normal(jac)
                    τ′ = _integration_measure(jac)
                    return K(xnode, (coords=y, normal=ν)) * L(ŷ) * τ′ * l′ * μ′
                end
                er ≤ max(atol) || @warn "hcubature failed to converge: (I,E) = $val, $er"
                return val
            end
            for (k, j) in enumerate(jglob)
                qx,qy = Xqnodes[i], Yqnodes[j]
                push!(out.I, i)
                push!(out.J, j)
                push!(out.V, W[k] - K(qx,qy)*weight(qy))
            end
        end
    end
    return out
end
