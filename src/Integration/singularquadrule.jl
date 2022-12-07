"""
    struct SingularQuadratureRule{D,Q,S} <: AbstractQuadratureRule{D}

A quadrature rule over `D` intended to integrate functions which are singular at
a known point `s ∈ D`.

A singular quadrature is rule is composed of a *regular* quadrature rule (e.g.
`GaussLegendre`) and a [`AbstractSingularityHandler`](@ref) to transform the regular
quadrature. The regular quadrature rule generates nodes and weights on the
`domain(sing_handler)`, and those are mapped into an appropriate quadrature over
`D = range(sing_handler)` using the singularity handler.

```jldoctest
import WavePropBase as WPB
q = WPB.Fejer(20)
shand = WPB.Kress()
qsing = WPB.SingularQuadratureRule(q, shand)

# check that we can accurately integrate the singular function `log(x)`
WPB.integrate(log,qsing) ≈ -1

# output

true
```

"""
struct SingularQuadratureRule{D,Q,S} <: AbstractQuadratureRule{D}
    qrule::Q
    singularity_handler::S
    function SingularQuadratureRule(q::Q, s::S) where {Q,S}
        @assert domain(s) == domain(q) "domain of quadrature must coincide with domain of singularity handler"
        D = typeof(image(s))
        return new{D,Q,S}(q, s)
    end
end

# getters
qrule(q::SingularQuadratureRule) = q.qrule
singularity_handler(q::SingularQuadratureRule) = q.singularity_handler

domain(qs::SingularQuadratureRule) = domain(singularity_handler(qs))
image(qs::SingularQuadratureRule) = image(singularity_handler(qs))

@generated function (qs::SingularQuadratureRule{D,Q,S})() where {D,Q,S}
    qstd = Q()
    shand = S()
    x̂, ŵ = qstd() # reference quadrature
    # use the change of variables cov to modify reference quadrature
    x = map(x -> shand(x), x̂)
    w = map(x̂, ŵ) do x̂, ŵ
        jac = jacobian(shand, x̂)
        μ = abs(det(jac))
        return μ * prod(ŵ)
    end
    return :($x, $w)
end

"""
    singular_quadrature(k,q::SingularQuadratureRule,s)

Return nodes and weights to integrate a function over `domain(q)` with a
factored weight `k`.
"""
function singular_quadrature(k, q::SingularQuadratureRule, s)
    x, w = q(s)
    T = Base.promote_op(k, eltype(x))
    assert_concrete_type(T)
    w = map(zip(x, w)) do (x, w)
        return k(x) * w
    end
    return x, w
end

"""
    singular_weights(k,xi,q::SingularQuadratureRule,s)
"""
function singular_weights(k, xi, q::SingularQuadratureRule, s)
    x, w = singular_quadrature(k, q, s)
    x = map(x -> x[1], x)
    ws = map(lag_basis) do li
        integrate(x, w) do x
            return f(x) * li(x)
        end
    end
    # wlag = barycentric_lagrange_weights(xi)
    # L    = barycentric_lagrange_matrix(xi,x,wlag)
    # return transpose(L)*w
    return ws
end

function singular_weights(k, qreg::AbstractQuadratureRule, qsin::SingularQuadratureRule, s)
    xq, wq = qsin(s)
    xi = qcoords(qreg)
    lag_basis = lagrange_basis(xi)
    map(lag_basis) do l
        integrate(xq, wq) do x
            return k(x) * l(x)
        end
    end
    # return transpose(L)*w
end

function singular_quadrature_suite(qreg::AbstractQuadratureRule{ReferenceLine},qsing::AbstractQuadratureRule{ReferenceLine})
    τ̂ = domain(qreg)
    @assert τ̂ == domain(qsing)
    xq = qcoords(qreg)
    nq = length(xq)
    nodes   = [SVector{1,Float64}[] for _ in 1:nq]
    weights = [Float64[] for _ in 1:nq]
    for (i,xi) in enumerate(xq)
        for el in decompose(τ̂,xi)
            xseg, wseg = qsing(el)
            append!(nodes[i], Vector(xseg))
            append!(weights[i], Vector(wseg))
        end
    end
    return nodes,weights
end

function integrate(f,qs::SingularQuadratureRule{<:Any,<:AdaptiveQuadrature,<:Any},args...;kwargs...)
    integrate_with_error(f,qs,args...;kwargs...)[1]
end

function integrate_with_error(f,qs::SingularQuadratureRule{<:Any,<:AdaptiveQuadrature,<:Any},args...)
    shand = singularity_handler(qs)
    q = qrule(qs)
    g = (x) -> f(shand(x))*integration_measure(shand,x)
    integrate_with_error(g,q,args...)
end
