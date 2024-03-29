"""
    struct PlotPoints

Structure used for dispatching `SVector` to plot recipes without type-piracy.
"""
struct PlotPoints end

@recipe function f(::PlotPoints, pts::SVector{N1,SVector{N2,Float64}}) where {N1,N2}
    if N2 == 2
        xx = SVector{N1,Float64}(pt[1] for pt in pts)
        yy = SVector{N1,Float64}(pt[2] for pt in pts)
        return xx, yy
    elseif N2 == 3
        xx = SVector{N1,Float64}(pt[1] for pt in pts)
        yy = SVector{N1,Float64}(pt[2] for pt in pts)
        zz = SVector{N1,Float64}(pt[3] for pt in pts)
        return xx, yy, zz
    else
        notimplemented()
    end
end

@recipe function f(::PlotPoints, pts::AbstractVector{<:SVector{N}}) where {N}
    if N == 2
        xx = [pt[1] for pt in pts]
        yy = [pt[2] for pt in pts]
        return xx, yy
    elseif N == 3
        xx = [pt[1] for pt in pts]
        yy = [pt[2] for pt in pts]
        zz = [pt[3] for pt in pts]
        return xx, yy, zz
    else
        notimplemented()
    end
end

@recipe function f(::PlotPoints, pts::AbstractMatrix{<:SVector})
    return PlotPoints(), vec(pts)
end

# plot a hyperrectangle
@recipe function f(rec::AbstractHyperRectangle{N}) where {N}
    seriestype := :path
    linecolor --> :black
    linestyle --> :solid
    label --> ""
    if N == 2
        pt1 = low_corner(rec)
        pt2 = high_corner(rec)
        x1, x2 = pt1[1], pt2[1]
        y1, y2 = pt1[2], pt2[2]
        @series SVector(x1, x2, x2, x1, x1), SVector(y1, y1, y2, y2, y1)
    elseif N == 3
        seriestype := :path
        pt1 = low_corner(rec)
        pt2 = high_corner(rec)
        x1, x2 = pt1[1], pt2[1]
        y1, y2 = pt1[2], pt2[2]
        z1, z2 = pt1[3], pt2[3]
        # upper and lower faces
        @series SVector(x1, x2, x2, x1, x1),
                SVector(y1, y1, y2, y2, y1),
                SVector(z1, z1, z1, z1, z1)
        @series SVector(x1, x2, x2, x1, x1),
                SVector(y1, y1, y2, y2, y1),
                SVector(z2, z2, z2, z2, z2)
        # lines connecting faces
        @series SVector(x1, x1), SVector(y1, y1), SVector(z1, z2)
        @series SVector(x2, x2), SVector(y1, y1), SVector(z1, z2)
        @series SVector(x2, x2), SVector(y2, y2), SVector(z1, z2)
        @series SVector(x1, x1), SVector(y2, y2), SVector(z1, z2)
    end
end

# plot domain --> plot all of its entities
@recipe function f(Ω::Domain)
    for ent in entities(Ω)
        @series begin
            ent
        end
    end
end

@recipe function f(Q::NystromMesh)
    @series begin
        aspect_ratio --> :equal
        mesh(Q)
    end
    @series begin
        label --> "Quadrature nodes"
        seriestype := :scatter
        aspect_ratio --> :equal
        PlotPoints(), collect(qcoords(Q))
    end
end

@recipe function f(mesh::AbstractMesh)
    label --> ""
    grid --> false
    aspect_ratio --> :equal
    for E in keys(mesh)
        for el in ElementIterator(mesh, E)
            @series begin
                el
            end
        end
    end
end

# FIXME: write better recipes
@recipe function f(el::LagrangeTriangle)
    label --> ""
    vtx = vals(el)
    for n in 1:3
        is = n
        ie = 1 + (n % 3)
        @series begin
            PlotPoints(), [vtx[is], vtx[ie]]
        end
    end
end

@recipe function f(el::LagrangeSquare)
    label --> ""
    vtx = vals(el)
    for n in 1:4
        is = n
        ie = 1 + (n % 4)
        @series begin
            PlotPoints(), [vtx[is], vtx[ie]]
        end
    end
end

@recipe function f(el::LagrangeLine)
    vtx = vals(el)
    return PlotPoints(), [vtx[1], vtx[2]]
end

"""
    struct PlotTree

Used to plot entire tree associated with a tree node, instead of just the node.
"""
struct PlotTree end

@recipe function f(::PlotTree, tree::ClusterTree; predicate=(x) -> isleaf(x))
    legend := false
    grid --> false
    aspect_ratio --> :equal
    # plot nodes
    blocks = filter_tree(predicate, tree)
    for block in blocks
        @series begin
            block
        end
    end
end

@recipe function f(tree::ClusterTree)
    legend := false
    grid --> false
    aspect_ratio --> :equal
    # plot points
    @series begin
        seriestype := :scatter
        markersize --> 2
        PlotPoints(), elements(tree)
    end
    # plot bounding box
    @series begin
        linestyle --> :solid
        seriescolor --> :black
        container(tree)
    end
end

@recipe function f(el::AbstractElement; npts=10)
    D = domain(el)
    grid --> false
    aspect_ratio --> :equal
    label --> ""
    if D isa ReferenceLine
        s = LinRange(0, 1, npts)
        pts = [el(SVector(v)) for v in s]
        x = [pt[1] for pt in pts]
        y = [pt[2] for pt in pts]
        x, y
    elseif D isa ReferenceSquare
        seriestype := :surface
        xrange = LinRange(0, 1, npts)
        yrange = LinRange(0, 1, npts)
        pts = [el(SVector(x, y)) for x in xrange, y in yrange]
        x = [pt[1] for pt in pts]
        y = [pt[2] for pt in pts]
        z = [pt[3] for pt in pts]
        x, y, z
    else
        notimplemented()
    end
end

# Plot recipes
@recipe function f(ent::ParametricEntity; step=0.1)
    D = domain(ent)
    grid --> false
    aspect_ratio --> :equal
    label --> ""
    if D isa HyperRectangle{1}
        msh = UniformCartesianMesh(D; step)
        iter = NodeIterator(msh)
        pts = [ent(v) for v in iter]
        x = [pt[1] for pt in pts]
        y = [pt[2] for pt in pts]
        x, y
    elseif D isa HyperRectangle{2}
        seriestype := :surface
        msh = UniformCartesianMesh(D; step)
        iter = NodeIterator(msh)
        pts = [ent(v) for v in iter]
        x = [pt[1] for pt in pts]
        y = [pt[2] for pt in pts]
        z = [pt[3] for pt in pts]
        x, y, z
    else
        notimplemented()
    end
end
