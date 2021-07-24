"""
    struct HyperRectangle{N,T}

Axis-aligned hyperrectangle in `N` dimensions given by `low_corner::SVector{N,T}` and
`high_corner::SVector{N,T}`.
"""
struct HyperRectangle{N,T}
    low_corner::SVector{N,T}
    high_corner::SVector{N,T}
end

low_corner(r::HyperRectangle)  = r.low_corner
high_corner(r::HyperRectangle) = r.high_corner

# promote
HyperRectangle(l::Tuple,h::Tuple)     = HyperRectangle(SVector(l),SVector(h))
HyperRectangle(l::SVector,h::SVector) = HyperRectangle(promote(l,h)...)
# 1d case
HyperRectangle(a::Number,b::Number) = HyperRectangle(SVector(a),SVector(b))

Base.:(==)(h1::HyperRectangle, h2::HyperRectangle) = (h1.low_corner == h2.low_corner) && (h1.high_corner == h2.high_corner)

Base.isapprox(h1::HyperRectangle,h2::HyperRectangle;kwargs...) = isapprox(h1.low_corner,h2.low_corner;kwargs...) && isapprox(h1.high_corner,h2.high_corner;kwargs...)

Base.in(point,h::HyperRectangle) = all(h.low_corner .<= point .<= h.high_corner)

Base.eltype(h::HyperRectangle{N,T}) where {N,T}     = T

ambient_dimension(h::HyperRectangle{N}) where {N}   = N

geometric_dimension(h::HyperRectangle{N}) where {N} = N

diameter(cub::HyperRectangle) = norm(high_corner(cub) .- low_corner(cub),2)

function bounding_box(pts,cube=false)
    isempty(pts)  && (error("data cannot be empty") )
    lb  = first(pts)
    ub  = first(pts)
    for pt in pts
        lb  = min.(lb,pt)
        ub  = max.(ub,pt)
    end
    if cube # fit a square/cube instead
        w  = maximum(ub-lb)
        xc = (ub + lb) / 2
        lb = xc .- w/2
        ub = xc .+ w/2
    end
    return HyperRectangle(lb,ub)
end
HyperRectangle(data,cube=false) = bounding_box(data,cube)

center(rec::HyperRectangle) = (rec.low_corner + rec.high_corner) / 2

radius(rec::HyperRectangle) = diameter(rec) / 2

"""
    split(rec::HyperRectangle,[axis]::Int,[place])

Split a hyperrectangle in two along the `axis` direction at the  position
`place`. Returns a tuple with the two resulting hyperrectangles.

When no `place` is given, defaults to splitting in the middle of the axis.

When no axis and no place is given, defaults to splitting along the largest axis.
"""
function Base.split(rec::HyperRectangle,axis,place)
    N            = ambient_dimension(rec)
    high_corner1 = svector(n-> n==axis ? place : rec.high_corner[n], N)
    low_corner2  = svector(n-> n==axis ? place : rec.low_corner[n], N)
    rec1         = HyperRectangle(rec.low_corner, high_corner1)
    rec2         = HyperRectangle(low_corner2,rec.high_corner)
    return (rec1, rec2)
end
function Base.split(rec::HyperRectangle,axis)
    place        = (rec.high_corner[axis] + rec.low_corner[axis])/2
    split(rec,axis,place)
end
function Base.split(rec::HyperRectangle)
    axis = argmax(rec.high_corner .- rec.low_corner)
    split(rec,axis)
end

"""
    distance(r1::HyperRectangle,r2::HyperRectangle)

The (minimal) Euclidean distance between a point `x ∈ r1` and `y ∈ r2`.
"""
function distance(rec1::HyperRectangle{N},rec2::HyperRectangle{N}) where {N}
    d2 = 0
    for i=1:N
        d2 += max(0,rec1.low_corner[i] - rec2.high_corner[i])^2 +
              max(0,rec2.low_corner[i] - rec1.high_corner[i])^2
    end
    return sqrt(d2)
end

"""
    distance(x::SVector,r::HyperRectangle)

The (minimal) Euclidean distance between the point `x` and any point `y ∈ r`.
"""
function distance(pt::SVector{N},rec::HyperRectangle{N}) where {N}
    d2 = 0
    for i=1:N
        d2 += max(0,pt[i] - rec.high_corner[i])^2 +
              max(0,rec.low_corner[i] - pt[i])^2
    end
    return sqrt(d2)
end
distance(rec::HyperRectangle{N},pt::SVector{N}) where {N} = distance(pt,rec)

function vertices(rec::HyperRectangle{2})
    lc = low_corner(rec)
    hc = high_corner(rec)
    SVector(
        SVector(lc[1],lc[2]),
        SVector(hc[1],lc[2]),
        SVector(hc[1],hc[2]),
        SVector(lc[1],hc[2])
    )
end
function vertices(rec::HyperRectangle{3})
    lc = low_corner(rec)
    hc = high_corner(rec)
    SVector(
        # lower face
        SVector(lc[1],lc[2],lc[3]),
        SVector(hc[1],lc[2],lc[3]),
        SVector(hc[1],hc[2],lc[3]),
        SVector(lc[1],hc[2],lc[3]),
        # upper face
        SVector(lc[1],lc[2],hc[3]),
        SVector(hc[1],lc[2],hc[3]),
        SVector(hc[1],hc[2],hc[3]),
        SVector(lc[1],hc[2],hc[3])
    )
end
