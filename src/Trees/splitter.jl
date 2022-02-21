"""
    abstract type AbstractSplitter

An `AbstractSplitter` is used to split a [`ClusterTree`](@ref). The interface
requires the following methods:
- `should_split(clt,splitter)` : return a `Bool` determining if the
  `ClusterTree` should be further divided
- `split!(clt,splitter)` : perform the splitting of the `ClusterTree` handling
  the necessary data sorting.

See [`GeometricSplitter`](@ref) for an example of an implementation.
"""
abstract type AbstractSplitter end

"""
    should_split(clt::ClusterTree,splitter::AbstractSplitter)

Determine whether or not a `ClusterTree` should be further divided.
"""
function should_split(clt, splitter::AbstractSplitter)
    abstract_method(splitter)
end

"""
    split!(clt::ClusterTree,splitter::AbstractSplitter)

Divide `clt` using the strategy implemented by `splitter`.
"""
function split!(clt, splitter::AbstractSplitter)
    abstract_method(splitter)
end

"""
    struct DyadicSplitter <: AbstractSplitter

Used to split an `N` dimensional `ClusterTree` into `2^N` children until at most
`nmax` points are contained in node.

## See also: [`AbstractSplitter`](@ref)
"""
Base.@kwdef struct DyadicSplitter <: AbstractSplitter
    nmax::Int = typemax(Int)
end

function should_split(node::ClusterTree, splitter::DyadicSplitter)
    return length(node) > splitter.nmax
end

function split!(parentcluster::ClusterTree, ::DyadicSplitter)
    d = ambient_dimension(parentcluster)
    clusters = [parentcluster]
    rec = container(parentcluster)
    rec_center = center(rec)
    for i = 1:d
        pos = rec_center[i]
        nel = length(clusters) #2^(i-1)
        for _ = 1:nel
            clt = popfirst!(clusters)
            append!(clusters, _binary_split!(clt, i, pos; parentcluster))
        end
    end
    return clusters
end

"""
    struct GeometricSplitter <: AbstractSplitter

Used to split a `ClusterTree` in half along the largest axis.
"""
Base.@kwdef struct GeometricSplitter <: AbstractSplitter
    nmax::Int = 50
end

should_split(node::ClusterTree, splitter::GeometricSplitter) = length(node) > splitter.nmax

function split!(cluster::ClusterTree, ::GeometricSplitter)
    rec = cluster.container
    wmax, imax = findmax(high_corner(rec) - low_corner(rec))
    left_node, right_node = _binary_split!(cluster, imax, low_corner(rec)[imax] + wmax / 2)
    return [left_node, right_node]
end

"""
    struct GeometricMinimalSplitter <: AbstractSplitter

Like [`GeometricSplitter`](@ref), but shrinks the children's containters.
"""
Base.@kwdef struct GeometricMinimalSplitter <: AbstractSplitter
    nmax::Int = 50
end

should_split(node::ClusterTree, splitter::GeometricMinimalSplitter) = length(node) > splitter.nmax

function split!(cluster::ClusterTree, ::GeometricMinimalSplitter)
    rec = cluster.container
    wmax, imax = findmax(high_corner(rec) - low_corner(rec))
    mid = low_corner(rec)[imax] + wmax / 2
    predicate = (x) -> x[imax] < mid
    left_node, right_node = _binary_split!(cluster, predicate)
    return [left_node, right_node]
end

"""
    struct PrincipalComponentSplitter <: AbstractSplitter
"""
Base.@kwdef struct PrincipalComponentSplitter <: AbstractSplitter
    nmax::Int = 50
end

should_split(node::ClusterTree, splitter::PrincipalComponentSplitter) = length(node) > splitter.nmax

function split!(cluster::ClusterTree, ::PrincipalComponentSplitter)
    pts = cluster._elements
    irange = cluster.index_range
    xc = center_of_mass(cluster)
    # compute covariance matrix for principal direction
    cov = sum(irange) do i
        x = coords(pts[i])
        (x - xc) * transpose(x - xc)
    end
    v = eigvecs(cov)[:, end]
    predicate = (x) -> dot(x - xc, v) < 0
    left_node, right_node = _binary_split!(cluster, predicate)
    return [left_node, right_node]
end

function center_of_mass(clt::ClusterTree)
    pts = clt._elements
    loc_idxs = clt.index_range
    # w    = clt.weights
    n = length(loc_idxs)
    # M    = isempty(w) ? n : sum(i->w[i],glob_idxs)
    # xc   = isempty(w) ? sum(i->pts[i]/M,glob_idxs) : sum(i->w[i]*pts[i]/M,glob_idxs)
    M = n
    xc = sum(i -> coords(pts[i]) / M, loc_idxs)
    return xc
end

"""
    struct CardinalitySplitter <: AbstractSplitter

Used to split a `ClusterTree` along the largest dimension if
`length(tree)>nmax`. The split is performed so the `data` is evenly distributed
amongst all children.

## See also: [`AbstractSplitter`](@ref)
"""
Base.@kwdef struct CardinalitySplitter <: AbstractSplitter
    nmax::Int = 50
end

should_split(node::ClusterTree, splitter::CardinalitySplitter) = length(node) > splitter.nmax

function split!(cluster::ClusterTree, ::CardinalitySplitter)
    points = cluster._elements
    irange = cluster.index_range
    rec = container(cluster)
    _, imax = findmax(high_corner(rec) - low_corner(rec))
    med = median(coords(points[i])[imax] for i in irange) # the median along largest axis `imax`
    predicate = (x) -> x[imax] < med
    left_node, right_node = _binary_split!(cluster, predicate)
    return [left_node, right_node]
end
