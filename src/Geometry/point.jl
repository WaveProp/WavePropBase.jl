# We will always reprenset a point as an `SVector`. Some aliases are defined
# below for convenience, but otherwise there is not much to defined regarding
# the representation of points.

"""
    const Point1D
    const Point1D(x1)

A point in 1D space, stored in a StaticArray.
Point1D = SVector{1, Float64}.
"""
const Point1D = SVector{1,Float64}

"""
    const Point2D
    const Point2D(x1, x2)
    const Point2D(x::NTuple{2, Float64})

A point in 2D space, stored in a StaticArray.
Point2D = SVector{2, Float64}.
"""
const Point2D = SVector{2,Float64}

"""
    const Point3D
    const Point3D(x1, x2, x3)
    const Point3D(x::NTuple{3, Float64})

A point in 3D space, stored in a StaticArray.
Point3D = SVector{3, Float64}.
"""
const Point3D = SVector{3,Float64}

"""
    const ComplexPoint3D
    const ComplexPoint3D(x1, x2, x3)
    const ComplexPoint3D(x::NTuple{3, ComplexF64})

A complex 3D point, stored in a StaticArray.
ComplexPoint3D = SVector{3, ComplexF64}.
"""
const ComplexPoint3D = SVector{3,ComplexF64}

"""
    const ComplexPoint2D
    const ComplexPoint2D(x1, x2)
    const ComplexPoint2D(x::NTuple{2, ComplexF64})

A complex 2D point, stored in a StaticArray.
ComplexPoint2D = SVector{2, ComplexF64}.
"""
const ComplexPoint2D = SVector{2,ComplexF64}
