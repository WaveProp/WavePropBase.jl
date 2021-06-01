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
        [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1]
    elseif N == 3
        seriestype := :path
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1], pt2[1]
        y1, y2 = pt1[2], pt2[2]
        z1, z2 = pt1[3], pt2[3]
        # upper and lower faces
        @series [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], [z1,z1,z1,z1,z1]
        @series [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], [z2,z2,z2,z2,z2]
        # lines connecting faces
        @series [x1,x1], [y1,y1], [z1,z2]
        @series [x2,x2], [y1,y1], [z1,z2]
        @series [x2,x2], [y2,y2], [z1,z2]
        @series [x1,x1], [y2,y2], [z1,z2]
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

# plot the mesh of a domain
@recipe function f(mesh::GenericMesh,Ω::Domain)
    view(mesh,Ω)
end
@recipe function f(mesh::SubMesh)
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

@recipe function f(el::LagrangeTriangle)
    label --> ""
    linecolor --> :black
    seriestype := :line
    vtx = vals(el)
    for n in 1:3
        is = n
        ie = 1 + (n%3)
        @series begin
            [vtx[is],vtx[ie]]
        end
    end
end

@recipe function f(el::LagrangeRectangle)
    label --> ""
    vtx = vals(el)
    for n in 1:4
        is = n
        ie = 1 + (n%4)
        @series begin
            [vtx[is],vtx[ie]]
        end
    end
end

@recipe function f(el::LagrangeLine)
    vtx = vals(el)
    [vtx[1],vtx[2]]
end

# plot a vector of points
@recipe function f(pts::AbstractVector{<:SVector{N}}) where {N}
    if N == 2
        xx = [pt[1] for pt in pts]
        yy = [pt[2] for pt in pts]
        return xx,yy
    elseif N == 3
        xx = [pt[1] for pt in pts]
        yy = [pt[2] for pt in pts]
        zz = [pt[3] for pt in pts]
        return xx,yy,zz
    else
        notimplemented()
    end
end
