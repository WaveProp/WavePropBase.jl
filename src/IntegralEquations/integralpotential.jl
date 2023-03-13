"""
    struct IntegralPotential

Represent a potential given by a `kernel` and a `source_mesh` over which
integration is performed.

`IntegralPotential`s are created using `IntegralPotential(kernel, source_mesh)`.

Evaluating an integral potential requires a density `σ` (defined over the
quadrature nodes of the source mesh) and a point `x` at which to evaluate the
integral
```math
\\int_{\\Gamma} K(\boldsymbol{x},\boldsymbol{y})\\sigma(y) ds_y, x \\not \\in \\Gamma
```

"""
struct IntegralPotential
    kernel::AbstractKernel
    source_mesh::AbstractMesh
end

kernel(pot::IntegralPotential) = pot.kernel
source_mesh(pot::IntegralPotential) = pot.source_mesh

function Base.getindex(pot::IntegralPotential, σ::AbstractVector)
    K = kernel(pot)
    Q = qnodes(source_mesh(pot))
    return (x) -> _evaluate_potential(K, σ, x, Q)
end

@noinline function _evaluate_potential(K, σ, x, Q)
    iter = zip(Q, σ)
    out = sum(iter) do (qi, σi)
        wi = weight(qi)
        return K(x, qi) * σi * wi
    end
    return out
end

function SingleLayerPotential(op::AbstractPDE, source)
    return IntegralPotential(SingleLayerKernel(op), source)
end
function DoubleLayerPotential(op::AbstractPDE, source)
    return IntegralPotential(DoubleLayerKernel(op), source)
end
