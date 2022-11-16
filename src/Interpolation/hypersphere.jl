struct HyperSphere{N,T}
    center::SVector{N,T}
    radius::T
end

center(s::HyperSphere) = s.center
radius(s::HyperSphere) = s.radius

# compute the circumscribed circle to a three dimensional triangle
# https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
function HyperSphere(el::AbstractElement{<:ReferenceTriangle})
    a = el((0, 0))
    b = el((1, 0))
    c = el((0, 1))
    xc = a +
         (norm(c - a)^2 * cross(cross(b - a, c - a), b - a) +
          norm(b - a)^2 * cross(cross(c - a, b - a), c - a)) /
         (2 * norm(cross(b - a, c - a))^2)
    r = norm(xc - a)
    return HyperSphere(xc, r)
end

# compute the circumcenter of a triangle in 3D
