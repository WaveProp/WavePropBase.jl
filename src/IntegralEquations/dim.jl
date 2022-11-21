Base.@kwdef struct DimParameters
    sources_oversample_factor::Float64 = 3
    sources_radius_multiplier::Float64 = 5
end

"""
    dim_correction(pde,X,Y,S,D[;location,p,derivative])

Given a `pde` and a (possibly innacurate) discretizations of its single and
double-layer operators `S` and `D` with domain `Y` and range `X`, compute
corrections `δS` and `δD` such that `S + δS` and `D + δD` are more accurate
approximations of the underlying integral operators.

The following optional keyword arguments are available:
- `location`: whether the target `X` lies on, inside, or outside the domain `Y`.
- `p::DimParameters`: parameters associated with the density interpolation
  method
- `derivative`: if true, compute the correction to the adjoint double-layer and
  hypersingular operators instead. In this case, `S` and `D` should contain a
  discretization of adjoint double-layer and hypersingular operators,
  respectively.
"""
function dim_correction(pde, X, Y, S, D; location=:onsurface, p=DimParameters(),
                        derivative=false)
    T = eltype(S)
    isvectorial = T <: SMatrix
    msg = "unrecognized value for kw `location`: received $location.
    Valid options are `:onsurface`, `:inside` and `:outside`."
    σ = location === :onsurface ? -0.5 :
        location === :inside ? -1 : location === :outside ? 0 : error(msg)
    dict_near = nearest_point_to_element(X, Y)
    # find first an appropriate set of source points to center the monopoles
    qmax = sum(size(etype2qtags(Y, E), 1) for E in keys(Y))
    ns = ceil(Int, p.sources_oversample_factor * qmax)
    pts = qcoords(Y)
    bbox = HyperRectangle(pts)
    xc = center(bbox)
    r = p.sources_radius_multiplier * radius(bbox)
    if ambient_dimension(Y) == 2
        xs = uniform_points_circle(ns, r, xc)
    elseif ambient_dimension(Y) == 3
        xs = fibonnaci_points_sphere(ns, r, xc)
    else
        notimplemented()
    end
    # compute traces of monopoles on the source mesh
    Xnodes = qnodes(X)
    Ynodes = qnodes(Y)
    G = SingleLayerKernel(pde, T)
    γ₁G = AdjointDoubleLayerKernel(pde, T)
    γ₀B = Matrix{T}(undef, length(Ynodes), ns)
    γ₁B = Matrix{T}(undef, length(Ynodes), ns)
    for k in 1:ns
        for j in 1:length(Ynodes)
            γ₀B[j, k] = G(Ynodes[j], xs[k])
            γ₁B[j, k] = γ₁G(Ynodes[j], xs[k])
        end
    end
    # integrate the monopoles/dipoles over Y with target on X. This is the
    # slowest step, and passing a custom S,D can accelerate this computation
    Θ = S * γ₁B - D * γ₀B
    for k in 1:ns
        for i in 1:length(Xnodes)
            if derivative
                Θ[i, k] += σ * γ₁G(Xnodes[i], xs[k])
            else
                Θ[i, k] += σ * G(Xnodes[i], xs[k])
            end
        end
    end
    # finally compute the corrected weights as sparse matrices
    Is, Js, Ss, Ds = Int[], Int[], T[], T[]
    dict_near = nearest_point_to_element(X, Y)
    for E in keys(Y)
        qtags = etype2qtags(Y, E)
        near_list = dict_near[E]
        nq, ne = size(qtags)
        @assert length(near_list) == ne
        M = Matrix{T}(undef, 2nq, ns)
        for n in 1:ne # loop over elements of type E
            # if there is nothing near, skip immediately to next element
            isempty(near_list[n]) && continue
            # the weights for points in near_list[n] must be corrected when
            # integrating over the current element
            jglob = @view qtags[:, n]
            M0 = @view γ₀B[jglob, :]
            M1 = @view γ₁B[jglob, :]
            copy!(view(M, 1:nq, :), M0)
            copy!(view(M, (nq + 1):(2nq), :), M1)
            F = qr!(blockmatrix_to_matrix(M))
            for i in near_list[n]
                Θi = @view Θ[i:i, :]
                W_ = (blockmatrix_to_matrix(Θi) / F.R) * adjoint(F.Q)
                W = matrix_to_blockmatrix(W_, T)
                for k in 1:nq
                    push!(Is, i)
                    push!(Js, jglob[k])
                    push!(Ss, -W[nq + k]) # single layer corresponds to α=0,β=-1
                    push!(Ds, W[k]) # double layer corresponds to α=1,β=0
                end
            end
        end
    end
    δS = sparse(Is, Js, Ss)
    δD = sparse(Is, Js, Ds)
    return δS, δD
end

"""
    nearest_point_to_element(X::NystroMesh,Y::NystromMesh; tol)

For each element `el`, return a list with the indices of all points in `X` for
which `el` is the nearest element. Ignore indices for which the distance exceeds `tol`.
"""
function nearest_point_to_element(X, Y::NystromMesh; tol=Inf)
    y = collect(qcoords(Y))
    kdtree = KDTree(y)
    dict = Dict(j => Int[] for j in 1:length(y))
    for i in eachindex(X)
        qtag, d = nn(kdtree, X[i])
        d > tol || push!(dict[qtag], i)
    end
    # dict[j] now contains indices in X for which the j quadrature node in Y is
    # the closest. Next we reverse the map
    etype2nearlist = Dict{DataType,Vector{Vector{Int}}}()
    for E in keys(Y)
        tags = etype2qtags(Y, E)
        nq, ne = size(tags)
        etype2nearlist[E] = nearlist = [Int[] for _ in 1:ne]
        for j in 1:ne # loop over each element of type E
            for q in 1:nq # loop over qnodes in the element
                qtag = tags[q, j]
                append!(nearlist[j], dict[qtag])
            end
        end
    end
    return etype2nearlist
end
function nearest_point_to_element(X::NystromMesh, Y::NystromMesh)
    if X === Y
        # when both surfaces are the same, the "near points" of an element are
        # simply its own quadrature points
        dict = Dict{DataType,Vector{Vector{Int}}}()
        for E in keys(Y)
            idx_dofs = etype2qtags(Y, E)
            dict[E] = map(i -> collect(i), eachcol(idx_dofs))
        end
    else
        dict = nearest_point_to_element(collect(qcoords(X)), Y)
    end
    return dict
end