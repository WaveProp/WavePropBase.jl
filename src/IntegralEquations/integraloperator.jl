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

function Base.Matrix(iop::IntegralOperator{T}) where {T}
    X = target_mesh(iop) |> qnodes
    Y = source_mesh(iop) |> qnodes
    K = kernel(iop)
    out = Matrix{T}(undef, length(X), length(Y))
    _iop_to_matrix!(out, K, X, Y)
    return out
end
@noinline function _iop_to_matrix!(out,K,X,Y)
    for i in 1:length(X)
        for j in 1:length(Y)
            out[i,j] = K(X[i], Y[j]) * weight(Y[j])
        end
    end
    return out
end

# convenience constructors
SingleLayerOperator(op::AbstractPDE, X, Y=X) = IntegralOperator(SingleLayerKernel(op), X, Y)
DoubleLayerOperator(op::AbstractPDE, X, Y=X) = IntegralOperator(DoubleLayerKernel(op), X, Y)
AdjointDoubleLayerOperator(op::AbstractPDE, X, Y=X) = IntegralOperator(AdjointDoubleLayerKernel(op), X, Y)
HyperSingularOperator(op::AbstractPDE, X, Y=X) = IntegralOperator(HyperSingularKernel(op), X, Y)


# Applying Laplace's double-layer to a constant will yield either 1 or -1,
# depending on whether the target point is inside or outside the obstacle
function isinside(x::SVector, mesh::NystromMesh, s=1)
    N   = ambient_dimension(mesh)
    pde = Laplace(dim=N)
    K   = DoubleLayerKernel(pde)
    u   = sum(qnodes(mesh)) do source
         K(x, source) * weight(source)
    end
    s*u + 0.5 < 0
    # u < 0
end
isinside(x::Tuple,mesh::NystromMesh) = isinside(SVector(x), mesh)
