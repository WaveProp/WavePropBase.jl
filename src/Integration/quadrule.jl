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
    qnodes(Y)

Return the quadrature nodes associated with `Y`.
"""
qnodes(q::AbstractQuadratureRule) = q()[1]

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
"""
struct Fejer{N} <: AbstractQuadratureRule{ReferenceLine} end

Fejer(n::Int) = Fejer{n}()

Fejer(; order::Int) = Fejer(order + 1)

# N point fejer quadrature integrates all polynomials up to degree N-1
order(::Fejer{N}) where {N} = N - 1

function (q::Fejer{N})() where {N}
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
    x, w = _get_gauss_qnodes_and_qweights(D, N)
    return :($x, $w)
end

"""
    TensorProductQuadrature{Q}

A tensor-product of one-dimension quadrature rules. Integrates over `[0,1]^d`,
where `d=length(quad)`.

# Examples
```julia
qx = Fejer(10)
qy = TrapezoidalOpen(15)
q  = TensorProductQuadrature(qx,qy)
```
"""
struct TensorProductQuadrature{Q} <: AbstractQuadratureRule{ReferenceSquare}
    quad::Q
end

# FIXME: this is a workaround the need to easily construct a tensor quadrature
# based only on the the types of the quadratures. Useful in generated functions,
# but there is probably a better way
function TensorProductQuadrature{Tuple{Q1,Q2}}() where {Q1,Q2}
    return TensorProductQuadrature(Q1(), Q2())
end

function TensorProductQuadrature(q...)
    return TensorProductQuadrature(q)
end

# FIXME: the current implementation is rather obscure. How should we handle the
# product quadrature rules in general? Also make this into a generated function.
function (q::TensorProductQuadrature)()
    N = length(q.quad)
    nodes = map(q -> q()[1], q.quad)
    weights = map(q -> q()[2], q.quad)
    nodes_iter = Iterators.product(nodes...)
    weights_iter = Iterators.product(weights...)
    x = map(x -> vcat(x...), nodes_iter)
    w = map(w -> prod(w), weights_iter)
    return SArray(x), w
end

"""
    qrule_for_reference_shape(ref,order)

Given a `ref`erence shape and a desired quadrature `order`, return
an appropiate quadrature rule.
"""
function qrule_for_reference_shape(ref, order)
    if ref isa ReferenceLine
        return Fejer(; order)
    elseif ref isa ReferenceSquare
        qx = qrule_for_reference_shape(ReferenceLine(), order)
        qy = qx
        return TensorProductQuadrature(qx, qy)
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
    qnodes::SVector{N,T}
    qweights::SVector{N,Float64}
end

function CustomQuadratureRule{D}(; qnodes, qweights) where {D}
    return CustomQuadratureRule{D,length(qnodes),eltype(qnodes)}(qnodes, qweights)
end
function CustomQuadratureRule(; domain::AbstractReferenceShape, qnodes, qweights)
    return CustomQuadratureRule{typeof(domain)}(; qnodes, qweights)
end

(q::CustomQuadratureRule)() = q.qnodes, q.qweights

const CustomLineQuadratureRule = CustomQuadratureRule{ReferenceLine}
const CustomTriangleQuadratureRule = CustomQuadratureRule{ReferenceTriangle}
const CustomSquareQuadratureRule = CustomQuadratureRule{ReferenceSquare}
const CustomTetrahedronQuadratureRule = CustomQuadratureRule{ReferenceTetrahedron}
