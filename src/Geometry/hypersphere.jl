"""
    struct HyperSphere{N,T}

A sphere in `ℜᴺ` described by its `center` and `radius`.
"""
struct HyperSphere{N,T}
    center::SVector{N,T}
    radius::T
end

center(s::HyperSphere) = s.center
radius(s::HyperSphere) = s.radius
