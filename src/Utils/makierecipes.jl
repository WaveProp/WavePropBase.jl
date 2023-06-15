import Makie as Makie

function tomakie_dim1(msh::AbstractMesh{N,T}) where {N,T}
    coords = Makie.Point{N,T}[]
    NAN    = svector(i->NaN*one(T),N)
    for E in keys(msh)
        iter = msh[E]
        D = domain(E)
        @assert D isa ReferenceLine
        vtxs = vertices(D)
        for el in iter
            for vtx in vtxs
                push!(coords, el(vtx))
            end
            # trick from Meshes.jl. Interleave with NaNs to plot segments
            push!(coords, NAN)
        end
    end
    return coords
end

function tomakie_dim2(msh::AbstractMesh{N,T}) where {N,T}
    coords = Makie.Point{N,T}[]
    connec = Int[]
    for E in keys(msh)
        iter = msh[E]
        D = domain(E)
        if D isa ReferenceTriangle
            vtxs = vertices(D)
            for el in iter
                for vtx in vtxs
                    push!(coords, el(vtx))
                    push!(connec, length(coords))
                end
            end
        elseif D isa ReferenceSquare
            # split square in two triangles for visualization
            vtxs_down = vertices(ReferenceTriangle())
            vtxs_up   = map(v-> -v .+ 1, vtx_down)
            for el in iter
                for vtxs in (vtxs_down, vtxs_up) # the two triangles
                    for vtx in vtxs
                        push!(coords, el(vtx))
                        push!(connec, length(coords))
                    end
                end
            end
        end
    end
    return coords,connec
end

function tomakie_dim3(msh::AbstractMesh{N,T}) where {N,T}
    coords = Makie.Point{N,T}[]
    connec = Int[]
    for E in keys(msh)
        iter = msh[E]
        D = domain(E)
        @assert D isa ReferenceTetrahedron
        vtxs = vertices(D)
        for el in iter
            for nf in 1:4 # four faces
                for (i,vtx) in enumerate(vtxs)
                    i == nf && continue # i-th face exclude the i-th vertex
                    push!(coords, el(vtx))
                    push!(connec, length(coords))
                end
            end
        end
    end
    return coords,connec
end

function Makie.convert_arguments(P::Type{<:Makie.Lines},msh::AbstractMesh)
    @assert geometric_dimension(msh) == 1 "Lines only supported for meshes of geometric dimension 1"
    coords = tomakie_dim1(msh)
    return (coords,)
end

function Makie.convert_arguments(P::Type{<:Makie.Poly},msh::AbstractMesh)
    gdim = geometric_dimension(msh)
    if gdim == 2
        coords, connec = tomakie_dim2(msh)
    elseif gdim == 3
        coords, connec = tomakie_dim3(msh)
    else
        error("Poly only supported for meshes of geometric dimension 2 or 3")
    end
    Makie.convert_arguments(P,coords,connec)
end

function Makie.convert_arguments(P::Type{<:Makie.Arrows},msh::AbstractMesh{N,T}) where {N,T}
    gdim = geometric_dimension(msh)
    adim = ambient_dimension(msh)
    codim = adim - gdim
    @assert codim == 1 "Arrows only supported for meshes of codimension 1"
    coords  = Makie.Point{N,T}[]
    normals = Makie.Point{N,T}[]
    for E in keys(msh)
        iter = msh[E]
        dom = domain(E)
        xc  = center(dom)
        for el in iter
            push!(coords, el(xc))
            push!(normals, normal(el,xc))
        end
    end
    Makie.convert_arguments(P,coords,normals)
end
