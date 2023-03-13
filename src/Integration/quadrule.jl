"""
    abstract type AbstractQuadratureRule{D}

A quadrature rule for integrating a function over the domain `D`.

Calling `q()` returns the nodes `x` and weights `w` for performing integration
over `domain(q)`.

Calling `q(el)` returns the nodes and weights for integrating over `el`.
"""
abstract type AbstractQuadratureRule{D} end

"""
    domain(f)

Given a function-like object `f: Ω → R`, return `Ω`.
"""
domain(q::AbstractQuadratureRule{D}) where {D} = D()

"""
    image(f)

Given a function-like object `f: Ω → R`, return `R`.
"""
function image(f)
    return abstractmethod(f)
end

"""
    qcoords(Y)

Return the coordinate of the quadrature nodes associated with `Y`.
"""
qcoords(q::AbstractQuadratureRule) = q()[1]

"""
    qweights(Y)

Return the quadrature weights associated with `Y`.
"""
qweights(q::AbstractQuadratureRule) = q()[2]

function (q::AbstractQuadratureRule)()
    return abstractmethod(q)
end

function (q::AbstractQuadratureRule)(el::AbstractElement)
    msg = "reference domain of quadrature rule and element must match"
    @assert domain(q) == domain(el) msg
    X̂, Ŵ = q()
    X = map(el, X̂)
    W = map(X̂, Ŵ) do x̂, ŵ
        jac = jacobian(el, x̂)
        μ = sqrt(det(jac' * jac))
        return μ * ŵ
    end
    return X, W
end

"""
    quadgen(el,q::AbstractQuadratureRule)

Return a set of points and weights for integrating over `el` based on the
quadrature rule `q`.
"""
quadgen(el::AbstractElement, q::AbstractQuadratureRule) = q(el)

"""
    integrate(f,q::AbstractQuadrature)
    integrate(f,x,w)

Integrate the function `f` using the quadrature rule `q`. This is simply
`sum(f.(x) .* w)`, where `x` and `w` are the quadrature nodes and weights, respectively.
"""
function integrate(f, q::AbstractQuadratureRule)
    x, w = q()
    if domain(q) == ReferenceLine()
        return integrate(x -> f(x[1]), x, w)
    else
        return integrate(f, x, w)
    end
end

function integrate(f, x, w)
    sum(zip(x, w)) do (x, w)
        return f(x) * prod(w)
    end
end

## Define some one-dimensional quadrature rules

"""
    struct Trapezoidal{N} <: AbstractQuadratureRule{ReferenceLine}

Closed `N`-point trapezoidal rule for integrating a function over the interval `[0,1]`.

# Examples:
```julia
q    = Trapezoidal(10)
f(x) = exp(x)*cos(x)
integrate(f,q)
```
"""
struct Trapezoidal{N} <: AbstractQuadratureRule{ReferenceLine} end

Trapezoidal(n::Int) = Trapezoidal{n}()

function (q::Trapezoidal{N})() where {N}
    h = 2 / (N - 1)
    x = [-1.0 + k * h for k in 0:(N - 1)]
    w = [h for _ in 1:N]
    w[1] = h / 2
    w[end] = h / 2
    # convert to static arrays
    xs = svector(i -> SVector(0.5 * (x[i] + 1)), N)
    ws = svector(i -> w[i] / 2, N)
    return xs, ws
end

"""
    struct TrapezoidalOpen{N} <: AbstractQuadratureRule{ReferenceLine}

Open trapezoidal rule.
"""
struct TrapezoidalOpen{N} <: AbstractQuadratureRule{ReferenceLine} end

TrapezoidalOpen(n::Int) = TrapezoidalOpen{n}()

# open trapezoidal rule
function _trapezoidal_open(n)
    h = 1 / n
    x = [(k - 0.5) * h for k in 1:n]
    w = [h for _ in 1:n]
    return x, w
end

function (q::TrapezoidalOpen{N})() where {N}
    x, w = _trapezoidal_open(N)
    # convert to static arrays
    xs = SVector{N}(SVector{1}.(x))
    ws = SVector{N}(w)
    return xs, ws
end

"""
    struct Fejer{N}

`N`-point Fejer's first quadrature rule for integrating a function over `[0,1]`.
Exactly integrates all polynomials of degree `≤ N-1`.

```jldoctest
import WavePropBase as WPB

q = WPB.Fejer(;order=10)

WPB.integrate(cos,q) ≈ sin(1) - sin(0)

# output

true
```

"""
struct Fejer{N} <: AbstractQuadratureRule{ReferenceLine} end

Fejer(n::Int) = Fejer{n}()

Fejer(; order::Int) = Fejer(order + 1)

# N point fejer quadrature integrates all polynomials up to degree N-1
order(::Fejer{N}) where {N} = N - 1

@generated function (q::Fejer{N})() where {N}
    theta = [(2j - 1) * π / (2 * N) for j in 1:N]
    x = -cos.(theta)
    w = zero(x)
    for j in 1:N
        tmp = 0.0
        for l in 1:floor(N / 2)
            tmp += 1 / (4 * l^2 - 1) * cos(2 * l * theta[j])
        end
        w[j] = 2 / N * (1 - 2 * tmp)
    end
    xs = svector(i -> SVector(0.5 * (x[i] + 1)), N)
    ws = svector(i -> w[i] / 2, N)
    return xs, ws
end

function integrate_fejer1d(f,
                           a=0,
                           b=1;
                           atol=0,
                           rtol=atol == 0 ? sqrt(eps(Float64)) : 0)
    w = b - a
    n = 5
    f̂ = (s) -> f(a + s * w) * w
    I1 = integrate(f̂, Fejer(n))
    I2 = integrate(f̂, Fejer(2 * n))
    er = norm(I2 - I1)
    while er > max(atol, I2 * rtol)
        n *= 2
        I1 = I2
        I2 = integrate(f̂, Fejer(2n))
        er = norm(I2 - I1)
        n > 100 && error("did not converge quickly enough")
    end
    return I2
end

"""
    struct GaussLegendre{N}

`N`-point Gauss-Legendre quadrature rule for integrating a function over `[0,1]`.
Exactly integrates all polynomials of degree `≤ 2N-1`.
"""
struct GaussLegendre{N} <: AbstractQuadratureRule{ReferenceLine} end

"""
    GaussLegendre(;order)

Construct a `GaussLegendre` of the desired order over the `[0,1]` interval.
"""
function GaussLegendre(; order=p)
    N = ceil(Int, (order + 1) / 2)
    return GaussLegendre{N}()
end

GaussLegendre(n::Int) = GaussLegendre{n}()

# N point Gauss quadrature integrates all polynomials up to degree 2N-1, yielding
# an error of order 2N
order(q::GaussLegendre{N}) where {N} = 2 * N - 1

# NOTE: the following code uses a Newton-Raphson method to find the roots of the
# Legendre polynomial of degree N. For low order, this is OK, but should not be
# used to high order. Better to rely on e.g. FastTranforms.jl for that matter.
@generated function (q::GaussLegendre{n})() where {n}
    # generate the nodes and weights on [-1,1]
    x = zeros(n)
    w = zeros(n)
    m = (n + 1) ÷ 2
    for i in 1:m
        z = cos(π * (i - 0.25) / (n + 0.5))
        z1 = z + 1
        pp = 0.0
        while abs(z - z1) > eps()
            p1 = 1.0
            p2 = 0.0
            for j in 1:n
                p3 = p2
                p2 = p1
                p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
            end
            pp = n * (z * p1 - p2) / (z^2 - 1)
            z1 = z
            z = z1 - p1 / pp
        end
        x[i] = -z
        x[n + 1 - i] = z
        w[i] = 2 / ((1 - z^2) * pp^2)
        w[n + 1 - i] = w[i]
    end
    # shift them to [0,1] and transform into a static vector
    xs = svector(i -> SVector(0.5 * (x[i] + 1)), n)
    ws = svector(i -> w[i] / 2, n)
    return xs, ws
end

"""
    struct Gauss{D,N} <: AbstractQuadratureRule{D}

Tabulated `N`-point symmetric Gauss quadrature rule for integration over `D`.
"""
struct Gauss{D,N} <: AbstractQuadratureRule{D}
    # gauss quadrature should be constructed using the order, and not the number
    # of nodes. This ensures you don't instantiate quadratures which are not
    # tabulated.
    function Gauss(; domain, order)
        domain == :triangle && (domain = ReferenceTriangle())
        domain == :tetrehedron && (domain = ReferenceTetrahedron())
        if domain isa ReferenceTriangle
            msg = "quadrature of order $order not available for ReferenceTriangle"
            haskey(TRIANGLE_GAUSS_ORDER_TO_NPTS, order) || error(msg)
            n = TRIANGLE_GAUSS_ORDER_TO_NPTS[order]
        elseif domain isa ReferenceTetrahedron
            msg = "quadrature of order $order not available for ReferenceTetrahedron"
            haskey(TETRAHEDRON_GAUSS_ORDER_TO_NPTS, order) || error(msg)
            n = TETRAHEDRON_GAUSS_ORDER_TO_NPTS[order]
        else
            error("Tabulated Gauss quadratures only available for `ReferenceTriangle` or `ReferenceTetrahedron`")
        end
        return new{typeof(domain),n}()
    end
end

function order(q::Gauss{ReferenceTriangle,N}) where {N}
    return TRIANGLE_GAUSS_NPTS_TO_ORDER[N]
end

function order(q::Gauss{ReferenceTetrahedron,N}) where {N}
    return TETRAHEDRON_GAUSS_NPTS_TO_ORDER[N]
end

@generated function (q::Gauss{D,N})() where {D,N}
    x, w = _get_gauss_qcoords_and_qweights(D, N)
    return :($x, $w)
end

"""
    TensorProductQuadrature{N,Q}

A tensor-product of one-dimension quadrature rules. Integrates over `[0,1]^N`.

# Examples
```julia
qx = Fejer(10)
qy = TrapezoidalOpen(15)
q  = TensorProductQuadrature(qx,qy)
```
"""
struct TensorProductQuadrature{N,Q} <: AbstractQuadratureRule{ReferenceHyperCube{N}}
    quads1d::Q
end

function TensorProductQuadrature(q...)
    N = length(q)
    Q = typeof(q)
    return TensorProductQuadrature{N,Q}(q)
end

function (q::TensorProductQuadrature{N})() where {N}
    nodes1d = ntuple(N) do i
        x1d, _ = q.quads1d[i]()
        return map(x -> x[1], x1d) # convert the `SVector{1,T}` to just `T`
    end
    weights1d = map(q -> q()[2], q.quads1d)
    nodes_iter = (SVector(x) for x in Iterators.product(nodes1d...))
    weights_iter = (prod(w) for w in Iterators.product(weights1d...))
    return nodes_iter, weights_iter
end

"""
    qrule_for_reference_shape(ref,order)

Given a `ref`erence shape and a desired quadrature `order`, return
an appropiate quadrature rule.
"""
function qrule_for_reference_shape(ref, order)
    if ref isa ReferenceLine
        return GaussLegendre(; order)
        # return Fejer(; order)
    elseif ref isa ReferenceSquare
        qx = qrule_for_reference_shape(ReferenceLine(), order)
        qy = qx
        return TensorProductQuadrature(qx, qy)
    elseif ref isa ReferenceCube
        qx = qrule_for_reference_shape(ReferenceLine(), order)
        qy = qz = qx
        return TensorProductQuadrature(qx, qy, qz)
    elseif ref isa ReferenceTriangle
        return Gauss(; domain=ref, order=order)
    elseif ref isa ReferenceTetrahedron
        return Gauss(; domain=ref, order=order)
    else
        error("no appropriate quadrature rule found.")
    end
end

"""
    struct CustomQuadratureRule{D,N,T} <: AbstractQuadratureRule{D}

`N`-point user-defined quadrature rule for integrating over `D`.
"""
struct CustomQuadratureRule{D<:AbstractReferenceShape,N,T<:SVector} <:
       AbstractQuadratureRule{D}
    qcoords::SVector{N,T}
    qweights::SVector{N,Float64}
end

function CustomQuadratureRule{D}(; qcoords, qweights) where {D}
    return CustomQuadratureRule{D,length(qcoords),eltype(qcoords)}(qcoords, qweights)
end
function CustomQuadratureRule(; domain::AbstractReferenceShape, qcoords, qweights)
    return CustomQuadratureRule{typeof(domain)}(; qcoords, qweights)
end

(q::CustomQuadratureRule)() = q.qcoords, q.qweights

const CustomLineQuadratureRule = CustomQuadratureRule{ReferenceLine}
const CustomTriangleQuadratureRule = CustomQuadratureRule{ReferenceTriangle}
const CustomSquareQuadratureRule = CustomQuadratureRule{ReferenceSquare}
const CustomTetrahedronQuadratureRule = CustomQuadratureRule{ReferenceTetrahedron}

# in the one-dimensional case
function lagrange_basis(q::AbstractQuadratureRule{ReferenceLine})
    nodes = qcoords(q)
    w = barycentric_lagrange_weights(q)
    N = length(w)
    # for some reason, I must pass a Val(N) to a function barrier so that the
    # anonymous function created is type stable and does not allocate
    return _lagrange_basis_1d(nodes, w, Val(N))
end

@noinline function _lagrange_basis_1d(nodes, w, ::Val{N}) where {N}
    (x) -> begin
        l = prod(xi -> x[1] - xi[1], nodes)
        svector(N) do j
            xj = nodes[j]
            num = l * w[j]
            den = x[1] - xj[1]
            return ifelse(den == 0, one(num), num / den)
        end
    end
end

function barycentric_lagrange_weights(q::AbstractQuadratureRule{ReferenceLine})
    # transform nodes to a vector
    nodes = Vector([x[1] for x in qcoords(q)])
    # comptue weights as a vector
    w = barycentric_lagrange_weights(nodes)
    ws = SVector{length(w)}(w)
    return ws
end

# in the two-dimensional case
function lagrange_basis(q::TensorProductQuadrature{2,<:Any})
    lx = lagrange_basis(q.quads1d[1])
    ly = lagrange_basis(q.quads1d[2])
    return _lagrange_basis2d(lx, ly)
end

@noinline function _lagrange_basis2d(lx, ly)
    return (x) -> begin
        vx = lx(SVector(x[1]))
        vy = ly(SVector(x[2]))
        vx * transpose(vy)
    end
end
