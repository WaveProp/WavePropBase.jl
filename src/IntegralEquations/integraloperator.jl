"""
    struct IntegralOperator{T} <: AbstractMatrix{T}

A discrete linear integral operator given by
```math
I[u](x) = \\int_{\\Gamma\\_s} K(x,y)u(y) ds_y, x \\in \\Gamma_{t}
```
where ``\\Gamma_s`` and ``\\Gamma_t`` are the source and target domains, respectively.
"""
struct IntegralOperator{T} <: AbstractMatrix{T}
    kernel::AbstractKernel
    target_mesh::AbstractMesh
    source_mesh::AbstractMesh
end

kernel(iop::IntegralOperator) = iop.kernel
target_mesh(iop::IntegralOperator) = iop.target_mesh
source_mesh(iop::IntegralOperator) = iop.source_mesh

function IntegralOperator(k, X, Y=X)
    T = return_type(k)
    msg = """IntegralOperator of nonbits being created"""
    isbitstype(T) || (@warn msg)
    return IntegralOperator{T}(k, X, Y)
end

function Base.size(iop::IntegralOperator)
    X = target_mesh(iop)
    Y = source_mesh(iop)
    return (length(qnodes(X)), length(qnodes(Y)))
end

function Base.getindex(iop::IntegralOperator, i::Integer, j::Integer)
    k = kernel(iop)
    targets = qnodes(target_mesh(iop))
    sources = qnodes(source_mesh(iop))
    return k(targets[i], sources[j]) * weight(sources[j])
end

# convenience constructors
SingleLayerOperator(op::AbstractPDE, X, Y=X) = IntegralOperator(SingleLayerKernel(op), X, Y)
DoubleLayerOperator(op::AbstractPDE, X, Y=X) = IntegralOperator(DoubleLayerKernel(op), X, Y)
function AdjointDoubleLayerOperator(op::AbstractPDE, X, Y=X)
    return IntegralOperator(AdjointDoubleLayerKernel(op), X, Y)
end
function HyperSingularOperator(op::AbstractPDE, X, Y=X)
    return IntegralOperator(HyperSingularKernel(op), X, Y)
end
