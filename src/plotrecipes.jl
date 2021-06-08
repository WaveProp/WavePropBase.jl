
"""
    struct PlotPoints

Structure used for dispatching plot recipes.
Its use also avoids the problem of type-piracy.
"""
struct PlotPoints end

# plot a vector of points
@recipe function f(::PlotPoints, pts::SVector{N1, SVector{N2, Float64}}) where {N1, N2}
    if N2 == 2
        xx = SVector{N1, Float64}(pt[1] for pt in pts)
        yy = SVector{N1, Float64}(pt[2] for pt in pts)
        return xx,yy
    elseif N2 == 3
        xx = SVector{N1, Float64}(pt[1] for pt in pts)
        yy = SVector{N1, Float64}(pt[2] for pt in pts)
        yy = SVector{N1, Float64}(pt[3] for pt in pts)
        return xx,yy,zz
    else
        notimplemented()
    end
end

# plot a hyperrectangle
@recipe function f(rec::HyperRectangle{N}) where {N}
    seriestype := :path
    linecolor --> :black
    linestyle --> :solid
    label --> ""
    if N == 2
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1], pt2[1]
        y1, y2 = pt1[2], pt2[2]
        @series SVector(x1,x2,x2,x1,x1), SVector(y1,y1,y2,y2,y1)
    elseif N == 3
        seriestype := :path
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1], pt2[1]
        y1, y2 = pt1[2], pt2[2]
        z1, z2 = pt1[3], pt2[3]
        # upper and lower faces
        @series SVector(x1,x2,x2,x1,x1), SVector(y1,y1,y2,y2,y1), SVector(z1,z1,z1,z1,z1)
        @series SVector(x1,x2,x2,x1,x1), SVector(y1,y1,y2,y2,y1), SVector(z2,z2,z2,z2,z2)
        # lines connecting faces
        @series SVector(x1,x1), SVector(y1,y1), SVector(z1,z2)
        @series SVector(x2,x2), SVector(y1,y1), SVector(z1,z2)
        @series SVector(x2,x2), SVector(y2,y2), SVector(z1,z2)
        @series SVector(x1,x1), SVector(y2,y2), SVector(z1,z2)
    end
end

# to plot domain a domain is to plot all of its entities. Assumes recipes for
# entities is available
@recipe function f(Ω::Domain)
    for ent in entities(Ω)
        @series begin
            ent
        end
    end
end

# plot all entities of the mesh
@recipe function f(mesh::GenericMesh)
    Ω = entities(mesh)|> Domain
    view(mesh,Ω)
end
# plot the mesh of a domain
@recipe function f(mesh::GenericMesh,Ω::Domain)
    view(mesh,Ω)
end
@recipe function f(mesh::SubMesh)
    label --> ""
    grid   --> false
    aspect_ratio --> :equal
    for E in keys(mesh)
        @series ElementIterator(mesh,E)
    end
end

@recipe function f(iter::ElementIterator)
    for el in iter
        @series begin
            el
        end
    end
end

# plot the mesh of a domain
@recipe function f(mesh::CartesianMesh)
    label --> ""
    grid   --> false
    aspect_ratio --> :equal
    for E in keys(mesh)
        for el in ElementIterator(mesh,E)
            @series begin
                el
            end
        end
    end
end

@recipe function f(el::LagrangeTriangle{3})
    label --> ""
    linecolor --> :black
    seriestype := :line
    n = vals(el)  # nodes
    @series PlotPoints(), SVector(n[1], n[2], n[3], n[1])
end
@recipe function f(el::LagrangeTriangle{6})
    label --> ""
    linecolor --> :black
    seriestype := :line
    n = vals(el)  # nodes
    # Reorder nodes for visualization purposes
    @series PlotPoints(), SVector(n[1], n[4], n[2], n[5], n[3], n[6], n[1])
end

@recipe function f(el::LagrangeRectangle)
    label --> ""
    n = vals(el)  # nodes
    @series PlotPoints(), SVector(n[1], n[2], n[3], n[4], n[1])
end

@recipe function f(el::LagrangeLine)
    n = vals(el)  # nodes
    @series PlotPoints(), SVector(n[1], n[2])
end


