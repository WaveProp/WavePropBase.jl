function vdim_correction(pde, X, Y, Γ, S, D, V; location, order)
    @assert X === Y
    N = ambient_dimension(pde)
    msg = "unrecognized value for kw `location`: received $location.
    Valid options are `:onsurface`, `:inside` and `:outside`."
    σ = location === :onsurface ? -0.5 :
        location === :inside ? -1.0 : location === :outside ? 0.0 : error(msg)
    @assert ambient_dimension(Y) == N "vdim only works for volume potentials"
    basis = _basis_vdim(pde, order)
    B, R = _vdim_auxiliary_quantities(basis, X, Y, Γ, σ, S, D, V)
    # compute sparse correction
    Is = Int[]
    Js = Int[]
    Vs = Float64[]
    nbasis = length(basis[1])
    max_cond = -Inf
    for E in keys(Y)
        qtags = etype2qtags(Y, E)
        nq, nel = size(qtags)
        @debug "nq,nbasis = $nq,$nbasis"
        for n in 1:nel
            # indices of nodes in element `n`
            j_idxs = qtags[:, n]
            i_idxs = j_idxs # FIXME: assumes that X === Y
            L = B[j_idxs, :] # vandermond matrix
            max_cond = max(max_cond, cond(L))
            for i in i_idxs
                wei = R[i:i, :] / L
                append!(Js, j_idxs)
                append!(Is, fill(i, length(j_idxs)))
                append!(Vs, wei)
            end
        end
    end
    @debug "maximum condition encountered: $max_cond"
    δV = sparse(Is, Js, Vs)
    return δV
end

function _vdim_auxiliary_quantities(basis, X, Y, Γ, σ, S, D, V)
    num_basis = length(basis[1])
    num_targets = length(qnodes(X))
    num_sources = length(qnodes(Y))
    num_bnd_sources = length(qnodes(Γ))
    # The compiler is not able to infer the type on the list comprehension
    # below, so we resort to preallocation and a regular loop
    # B   = [b(x) for x in dofs(Y), b in basis[1]]
    # γ₀B = [b(x) for x in dofs(Γ), b in basis[2]]
    # γ₁B = [b(x) for x in dofs(Γ), b in basis[3]]
    B = Matrix{Float64}(undef, num_sources, num_basis)
    γ₀B = Matrix{Float64}(undef, num_bnd_sources, num_basis)
    γ₁B = Matrix{Float64}(undef, num_bnd_sources, num_basis)
    for n in 1:num_basis
        for i in 1:num_sources
            y = qnodes(Y)[i]
            B[i, n] = basis[1][n](y)
        end
        for i in 1:num_bnd_sources
            y = qnodes(Γ)[i]
            γ₀B[i, n] = basis[2][n](y)
            γ₁B[i, n] = basis[3][n](y)
        end
    end
    R = S * γ₁B - D * γ₀B - V * B
    for i in 1:num_targets
        for n in 1:num_basis
            x = coords(qnodes(X)[i])
            R[i, n] += σ * basis[2][n](x)
        end
    end
    return B, R
end

"""
    _basis_vdim(ℒ::AbstractPDE,ℙₖ::PolynomialSpace)

Return a set of polynomial `pₙ` such that `ℒ[pₙ]` gives the monomial basis for
`ℙₖ`.
"""
function _basis_vdim(::Laplace{2}, order)
    # define the monomial basis functions, and the corresponding solutions
    # P0
    r_00 = (dof) -> 1.0
    p_00 = (dof) -> begin
        x = coords(dof)
        1 / 4 * (x[1]^2 + x[2]^2)
    end
    dp_00 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 2 * dot(x, n)
    end
    R, P, dP = (r_00,), (p_00,), (dp_00,)
    order == 0 && return R, P, dP
    # P1
    r_10 = (dof) -> coords(dof)[1]
    p_10 = (dof) -> begin
        x, y = coords(dof)
        1 / 6 * x^3
    end
    dp_10 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 2 * x[1]^2 * n[1]
    end
    r_01 = (dof) -> coords(dof)[2]
    p_01 = (dof) -> begin
        x = coords(dof)
        1 / 6 * x[2]^3
    end
    dp_01 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 2 * x[2]^2 * n[2]
    end
    R, P, dP = (R..., r_10, r_01), (P..., p_10, p_01), (dP..., dp_10, dp_01)
    order == 1 && return R, P, dP
    # P2
    r_20 = (dof) -> coords(dof)[1]^2
    p_20 = (dof) -> begin
        x, y = coords(dof)
        1 / 12 * x^4
    end
    dp_20 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 3 * x[1]^3 * n[1]
    end
    r_11 = (dof) -> coords(dof)[1] * coords(dof)[2]
    p_11 = (dof) -> begin
        x, y = coords(dof)
        x * y / 12 * (x^2 + y^2)
    end
    dp_11 = (dof) -> begin
        x, y = coords(dof)
        nx, ny = normal(dof)
        ((3 * x^2 * y + y^3) * nx + (x^3 + 3 * x * y^2) * ny) / 12
    end
    r_02 = (dof) -> coords(dof)[2]^2
    p_02 = (dof) -> begin
        x, y = coords(dof)
        1 / 12 * y^4
    end
    dp_02 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 3 * x[2]^3 * n[2]
    end
    R, P, dP = (R..., r_20, r_11, r_02), (P..., p_20, p_11, p_02),
               (dP..., dp_20, dp_11, dp_02)
    order == 2 && return R, P, dP
    # the rest has not been implemented
    return notimplemented()
end

function _basis_vdim(::Laplace{3}, order)
    # define the monomial basis functions, and the corresponding solutions
    # P0
    r_000 = (dof) -> 1.0
    p_000 = (dof) -> begin
        x = coords(dof)
        1 / 6 * (x[1]^2 + x[2]^2 + x[3]^2)
    end
    dp_000 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 3 * dot(x, n)
    end
    R, P, dP = (r_000,), (p_000,), (dp_000,)
    order == 0 && return R, P, dP
    # P1
    r_100 = (dof) -> coords(dof)[1]
    p_100 = (dof) -> begin
        x = coords(dof)
        1 / 6 * x[1]^3
    end
    dp_100 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 2 * x[1]^2 * n[1]
    end
    r_010 = (dof) -> coords(dof)[2]
    p_010 = (dof) -> begin
        x = coords(dof)
        1 / 6 * x[2]^3
    end
    dp_010 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 2 * x[2]^2 * n[2]
    end
    r_001 = (dof) -> coords(dof)[3]
    p_001 = (dof) -> begin
        x = coords(dof)
        1 / 6 * x[3]^3
    end
    dp_001 = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        1 / 2 * x[3]^2 * n[3]
    end
    R, P, dP = (R..., r_100, r_010, r_001), (P..., p_100, p_010, p_001),
               (dP..., dp_100, dp_010, dp_001)
    order == 1 && return R, P, dP
    return notimplemented()
    # P2
end
