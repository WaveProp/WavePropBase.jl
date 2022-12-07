Base.@kwdef struct AdaptiveQuadrature{D,Q<:EmbeddedQuadratureRule} <:
                   AbstractQuadratureRule{D}
    qrule::Q = Q()
    atol::Float64 = 0
    rtol::Float64 = atol > 0 ? 0 : sqrt(eps())
    maxsplit::Int = 100 # reasonable?
end

function (::AdaptiveQuadrature)()
    return error("calling `q()` for an adaptive quadrature `q` is not allowed. Use `integrate(f,q)` instead.")
end

function integrate(f, q::AdaptiveQuadrature, args...; kwargs...)
    return integrate_with_error(f, q, args...; kwargs...)[1]
end

function integrate_with_error(f::F, q::AdaptiveQuadrature,
                              heap=allocate_buffer(f, q)) where {F}
    I, E = integrate_with_error(f, q.qrule)
    nsplit = 0
    # a quick check to see if splitting is really needed
    if E < q.atol || E < q.rtol * norm(I) || nsplit >= q.maxsplit
        return I, E
    end
    # split is needed, so push the element to the heap and begin
    empty!(heap.valtree)
    T = eltype(heap).parameters[1]
    t = T(vertices(domain(q)))
    push!(heap, t => (I, E))
    I, E = _integrate_with_error!(f, heap, I, E, nsplit, q)
    return I, E
end
function _integrate_with_error!(f::F, heap, I, E, nsplit, q) where {F}
    while E > q.atol && E > q.rtol * norm(I) && nsplit < q.maxsplit
        t, (Ic, Ec) = pop!(heap)
        I -= Ic
        E -= Ec
        for child in split(t, q)
            # since the jacobian is constant, factor it out of the integration
            μ = integration_measure(child)
            g = (s) -> f(child(s))
            v = integrate_with_error(g, q.qrule)
            Inew, Enew = μ .* v
            I += Inew
            E += Enew
            push!(heap, child => (Inew, Enew))
        end
        nsplit += 1
    end
    nsplit >= q.maxsplit && @warn "maximum number of steps reached: $I, $E"
    return I, E
end

function allocate_buffer(f, q::AdaptiveQuadrature)
    # try to infer the type of the element that will be returned by f
    x, w = q.qrule()
    T = eltype(w) #
    S = return_type(f, eltype(x))
    TS = Base.promote_op(*, S, T)
    isbitstype(TS) || (@warn "non bitstype detected")
    d = domain(q)
    if d isa ReferenceLine
        D = Line1D{T}
    elseif d isa ReferenceTriangle
        D = Triangle2D{T}
    else
        error("unable to allocate buffer: domain of quadrature rule not supported")
    end
    # the heap of adaptive quadratures have elemtns of the form el => (I,E),
    # where I and E are the value and error estimate over the element el. The
    # ordering used is based on the negative of the error
    ord = Base.Order.By(el -> -el[2][2])
    heap = BinaryHeap{Pair{D,Tuple{TS,T}}}(ord)
    return heap
end


# some specific adaptive quadratures
const CubTri = AdaptiveQuadrature{ReferenceTriangle,
                                  EmbeddedQuadratureRule{Radon5,Laurie8,ReferenceTriangle}}

function Base.split(t::Triangle2D{T}, q::CubTri) where {T}
    p1, p2, p3 = vals(t)
    m1 = (p1 + p2) / 2
    m2 = (p2 + p3) / 2
    m3 = (p1 + p3) / 2
    t1 = Triangle2D{T}(p1, m1, m3)
    t2 = Triangle2D{T}(m1, p2, m2)
    t3 = Triangle2D{T}(m2, p3, m3)
    t4 = Triangle2D{T}(m1, m2, m3)
    return t1, t2, t3, t4
end

const G7K15 = AdaptiveQuadrature{ReferenceLine,
                                 EmbeddedQuadratureRule{GaussLegendre{7},Kronrod{15},
                                                        ReferenceLine}}

function Base.split(l::Line1D{T}, ::G7K15) where {T}
    p1, p2 = vals(l)
    m = (p1 + p2) / 2
    l1 = Line1D{T}(p1, m)
    l2 = Line1D{T}(m, p2)
    return l1, l2
end
