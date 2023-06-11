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
    hcubature_correction(iop::IntegralOperator; maxdist, kwargs...)

Given an integral operator `iop`, this function provides a sparse correction to
`iop` for the entries `i,j` such that the distance between the `i`-th target and
the `j`-th source is less than `maxdist`.

The correction is computed by using an adaptive quadrature based on the
`hcubature` routine from `HCubature.jl`; the `kwargs` are passed to `hcubature`,
and (see [its documentation](https://github.com/JuliaMath/HCubature.jl) for more
details).
"""
function hcubature_correction(iop::IntegralOperator; maxdist, atol = 1e-8, maxevals = 1000, initdiv = 1)
    # unpack type-unstable fields in iop, allocate output, and dispatch
    X, Y, K = target_mesh(iop), source_mesh(iop), kernel(iop)
    dict_near = near_interaction_list(X, Y; tol=maxdist)
    T  = eltype(iop)
    out = (I=Int[], J=Int[], V=T[])
    for E in keys(Y)
        # dispatch on element type
        nearlist = dict_near[E]
        qreg = etype2qrule(Y, E)
        L = lagrange_basis(qreg)
        iter = Y[E]
        _hcubature_correction_etype!(out, iter, qreg, L, nearlist, X, Y, K, maxdist, atol, maxevals, initdiv)
    end
    m,n = size(iop)
    return sparse(out.I,out.J,out.V, m, n)
end

@noinline function _hcubature_correction_etype!(out,iter,qreg,L,nearlist,X,Y,K,maxdist, atol, maxevals, initdiv)
    E = eltype(iter)
    Xqnodes = qnodes(X)
    Yqnodes = qnodes(Y)
    τ̂ = domain(E)
    N = geometric_dimension(τ̂)
    a, b = ntuple(i -> 0, N), ntuple(i -> 1, N)
    τ̂ === ReferenceLine() || notimplemented()
    x̂, ŵ = collect.(qreg())
    el2qtags = etype2qtags(Y, E)
    buffer = hcubature_buffer(x -> one(eltype(out.V)) * L(a) * first(ŵ), a, b)
    success = true
    maxer   = 0.0
    for (n, el) in enumerate(iter)
        jglob = view(el2qtags, :, n)
        inear = nearlist[n]
        for i in inear
            xnode = Xqnodes[i]
            # closest quadrature node
            dmin, j = findmin(n -> norm(coords(xnode) - coords(Yqnodes[jglob[n]])),1:length(jglob))
            x̂nearest = x̂[j]
            dmin > maxdist && continue
            issingular = iszero(dmin)
            # use hcubature for singular integration of lagrange basis
            if issingular
                W1, er1 = hcubature(a, Tuple(x̂nearest); buffer, atol = atol/2, maxevals, initdiv) do ŷ
                    y  = el(ŷ)
                    jac = jacobian(el, ŷ)
                    ν = _normal(jac)
                    τ′ = _integration_measure(jac)
                    return K(xnode, (coords=y, normal=ν)) * L(ŷ) * τ′
                end
                W2, er2 = hcubature(Tuple(x̂nearest),b; buffer, atol = atol/2, maxevals, initdiv) do ŷ
                    y  = el(ŷ)
                    jac = jacobian(el, ŷ)
                    ν = _normal(jac)
                    τ′ = _integration_measure(jac)
                    return K(xnode, (coords=y, normal=ν)) * L(ŷ) * τ′
                end
                W  = W1 + W2
                er = er1 + er2
            else
                W, er = hcubature(a, b; buffer, atol, maxevals, initdiv) do ŷ
                    y  = el(ŷ)
                    y == coords(xnode) && (@warn "possibly singular point")
                    jac = jacobian(el, ŷ)
                    ν = _normal(jac)
                    τ′ = _integration_measure(jac)
                    return K(xnode, (coords=y, normal=ν)) * L(ŷ) * τ′
                end
            end
            er ≤ max(atol) || (success = false)
            maxer = max(maxer, er)
            for (k, j) in enumerate(jglob)
                qx,qy = Xqnodes[i], Yqnodes[j]
                push!(out.I, i)
                push!(out.J, j)
                push!(out.V, W[k] - K(qx,qy)*weight(qy))
            end
        end
    end
    msg = """hcubature failed to converge to tolerance $atol within $maxevals
    function evaluations: maximum error was $maxer. Consider increasing the
    `atol` and/or `maxevals`"""
    success || @warn msg
    return out
end
