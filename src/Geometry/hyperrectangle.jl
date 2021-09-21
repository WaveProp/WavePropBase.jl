
"""
    abstract type AbstractHyperRectangle{N,T}

Axis-aligned hyperrectangle in `N` dimensions with coordinates
of type `SVector{N,T}`.
"""
abstract type AbstractHyperRectangle{N,T} end

Base.eltype(::AbstractHyperRectangle{N,T}) where {N,T}     = T
ambient_dimension(::AbstractHyperRectangle{N}) where {N}   = N
geometric_dimension(::AbstractHyperRectangle{N}) where {N} = N
low_corner(r::AbstractHyperRectangle)  = abstractmethod(r)
high_corner(r::AbstractHyperRectangle) = abstractmethod(r)
diameter(r::AbstractHyperRectangle)    = abstractmethod(r)
center(r::AbstractHyperRectangle)      = abstractmethod(r)
radius(r::AbstractHyperRectangle)      = abstractmethod(r)
vertices(r::AbstractHyperRectangle)    = abstractmethod(r)

Base.in(point,h::AbstractHyperRectangle) = all(low_corner(h) .<= point .<= high_corner(h))

function vertices(rec::AbstractHyperRectangle{2})
    lc = low_corner(rec)
    hc = high_corner(rec)
    return SVector(
        SVector(lc[1],lc[2]),
        SVector(hc[1],lc[2]),
        SVector(hc[1],hc[2]),
        SVector(lc[1],hc[2])
    )
end
function vertices(rec::AbstractHyperRectangle{3})
    lc = low_corner(rec)
    hc = high_corner(rec)
    return SVector(
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

"""
    distance(r1::AbstractHyperRectangle,r2::AbstractHyperRectangle)

The (minimal) Euclidean distance between a point `x ∈ r1` and `y ∈ r2`.
"""
function distance(rec1::AbstractHyperRectangle{N},rec2::AbstractHyperRectangle{N}) where {N}
    d2 = 0
    rec1_low_corner  = low_corner(rec1)
    rec1_high_corner = high_corner(rec1)
    rec2_low_corner  = low_corner(rec2)
    rec2_high_corner = high_corner(rec2)
    for i=1:N
        d2 += max(0,rec1_low_corner[i] - rec2_high_corner[i])^2 +
              max(0,rec2_low_corner[i] - rec1_high_corner[i])^2
    end
    return sqrt(d2)
end

"""
    distance(x::SVector,r::HyperRectangle)

The (minimal) Euclidean distance between the point `x` and any point `y ∈ r`.
"""
function distance(pt::SVector{N},rec::AbstractHyperRectangle{N}) where {N}
    d2 = 0
    rec_low_corner  = low_corner(rec)
    rec_high_corner = high_corner(rec)
    for i=1:N
        d2 += max(0,pt[i] - rec_high_corner[i])^2 +
              max(0,rec_low_corner[i] - pt[i])^2
    end
    return sqrt(d2)
end
distance(rec::AbstractHyperRectangle{N},pt::SVector{N}) where {N} = distance(pt,rec)

######
# HyperRectangle
######

"""
    struct HyperRectangle{N,T}

Axis-aligned hyperrectangle in `N` dimensions given by `low_corner::SVector{N,T}` and
`high_corner::SVector{N,T}`.
"""
struct HyperRectangle{N,T} <: AbstractHyperRectangle{N,T}
    low_corner::SVector{N,T}
    high_corner::SVector{N,T}
end
HyperRectangle(l::Tuple,h::Tuple)     = HyperRectangle(SVector(l),SVector(h))
HyperRectangle(l::SVector,h::SVector) = HyperRectangle(promote(l,h)...)
# 1d case
HyperRectangle(a::Number,b::Number) = HyperRectangle(SVector(a),SVector(b))

low_corner(r::HyperRectangle)  = r.low_corner
high_corner(r::HyperRectangle) = r.high_corner
center(r::HyperRectangle)   = (low_corner(r) + high_corner(r)) / 2
diameter(r::HyperRectangle) = norm(high_corner(r) .- low_corner(r),2)
radius(r::HyperRectangle)   = diameter(r) / 2
Base.isapprox(h1::HyperRectangle,h2::HyperRectangle;kwargs...) = isapprox(h1.low_corner,h2.low_corner;kwargs...) && isapprox(h1.high_corner,h2.high_corner;kwargs...)

######
# HyperCube
######

"""
    struct HyperCube{N,T}

Axis-aligned hypercube in `N` dimensions given by `low_corner::SVector{N,T}` and
`side::T`.
"""
struct HyperCube{N,T} <: AbstractHyperRectangle{N,T}
    low_corner::SVector{N,T}
    side::T
end
# 1d case
HyperCube(low_corner::Number,side::Number) = HyperCube(SVector(low_corner),side)

side(r::HyperCube) = r.side
low_corner(r::HyperCube)  = r.low_corner
high_corner(r::HyperCube) = low_corner(r) .+ side(r)
center(r::HyperCube) = low_corner(r) .+ side(r)/2
diameter(r::HyperCube{N}) where N = side(r)*sqrt(N)
radius(r::HyperCube) = diameter(r) / 2
Base.isapprox(h1::HyperCube,h2::HyperCube;kwargs...) = isapprox(h1.low_corner,h2.low_corner;kwargs...) && isapprox(h1.side,h2.side;kwargs...)

######
# Utils
######

function bounding_box(els,cube=false)
    isempty(els)  && (error("data cannot be empty") )
    lb  = first(els) |> coords
    ub  = first(els) |> coords
    for el in els
        pt = coords(el)
        lb  = min.(lb,pt)
        ub  = max.(ub,pt)
    end
    if cube # fit a square/cube instead
        w  = maximum(ub-lb)
        xc = (ub + lb) / 2
        lb = xc .- w/2
        ub = xc .+ w/2
        # TODO: return HyperCube instead
    end
    return HyperRectangle(lb,ub)
end
HyperRectangle(data,cube=false) = bounding_box(data,cube)

"""
    split(rec::AbstractHyperRectangle,[axis]::Int,[place])

Split a hyperrectangle in two along the `axis` direction at the  position
`place`. Returns a tuple with the two resulting hyperrectangles.

When no `place` is given, defaults to splitting in the middle of the axis.

When no axis and no place is given, defaults to splitting along the largest axis.
"""
function Base.split(rec::AbstractHyperRectangle,axis,place)
    rec_low_corner  = low_corner(rec)
    rec_high_corner = high_corner(rec)
    N            = ambient_dimension(rec)
    high_corner1 = svector(n-> n==axis ? place : rec_high_corner[n], N)
    low_corner2  = svector(n-> n==axis ? place : rec_low_corner[n], N)
    rec1         = HyperRectangle(rec_low_corner, high_corner1)
    rec2         = HyperRectangle(low_corner2,rec_high_corner)
    return (rec1, rec2)
end
function Base.split(rec::AbstractHyperRectangle,axis)
    place = (high_corner(rec)[axis] + low_corner(rec)[axis])/2
    return split(rec,axis,place)
end
function Base.split(rec::AbstractHyperRectangle)
    axis = argmax(high_corner(rec) .- low_corner(rec))
    return split(rec,axis)
end