# FIXME: document this
struct IntegralPotential
    kernel::AbstractKernel
    source_mesh::AbstractMesh
end

kernel(pot::IntegralPotential) = pot.kernel
source_mesh(pot::IntegralPotential) = pot.source_mesh

function (pot::IntegralPotential)(σ::AbstractVector, x)
    f = kernel(pot)
    Γ = source_mesh(pot)
    iter = zip(qnodes(Γ), σ)
    out = sum(iter) do (source, σ)
        w = weight(source)
        return f(x, source) * σ * w
    end
    return out
end

Base.getindex(pot::IntegralPotential, σ::AbstractVector) = (x) -> pot(σ, x)

function SingleLayerPotential(op::AbstractPDE, source)
    return IntegralPotential(SingleLayerKernel(op), source)
end
function DoubleLayerPotential(op::AbstractPDE, source)
    return IntegralPotential(DoubleLayerKernel(op), source)
end
