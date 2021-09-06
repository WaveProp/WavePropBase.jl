"""
    mutable struct ClusterTree{T,S,D}

Tree structure used to cluster elements of type `T` into containers of type `S`.
An additional `data` field of type `D` can be associated with each node (it
defaults to `D=Nothing`).

# Fields:
- `_elements::T` : vector containing the sorted elements.
- `loc_idxs::UnitRange{Int}` : indices of elements contained in the current node.
- `container::S` : container for the elements in the current node.
- `loc2glob::Vector{Int}` : permutation from global local indexing system to the
  original (global) indexing system used before the construction of the tree.
- `children::Vector{ClusterTree{N,T,D}}`
- `parent::ClusterTree{N,T,D}`
- `data::D` : generic data field parametrically typed on `D`.
"""
mutable struct ClusterTree{T,S,D}
    _elements::Vector{T}
    container::S
    index_range::UnitRange{Int}
    loc2glob::Vector{Int}
    children::Vector{ClusterTree{T,S,D}}
    parent::ClusterTree{T,S,D}
    data::D
    # inner constructors handling missing fields
    function ClusterTree{D}(els::Vector{T},container::S,loc_idxs,loc2glob,children,parent,data=nothing) where {T,S,D}
        clt = new{T,S,D}(els,container,loc_idxs,loc2glob)
        clt.children = isnothing(children) ? Vector{typeof(clt)}() : children
        clt.parent   = isnothing(parent)   ? clt : parent
        clt.data     = isnothing(data)     ? D() : data
        return clt
    end
end

# convenience functions and getters
root_elements(clt::ClusterTree) = clt._elements
index_range(clt::ClusterTree) = clt.index_range
children(clt::ClusterTree) = clt.children
parent(clt::ClusterTree)   = clt.parent
container(clt::ClusterTree) = clt.container
elements(clt::ClusterTree) = view(root_elements(clt),clt.index_range)
loc2glob(clt::ClusterTree) = clt.loc2glob

containter_type(clt::ClusterTree{T,S}) where {T,S} = S
element_type(clt::ClusterTree{T}) where {T} = T

# interface to AbstractTrees. No children is determined by an empty tuple for
# AbstractTrees.
AbstractTrees.children(t::ClusterTree) = isleaf(t) ? () : t.children
AbstractTrees.nodetype(t::ClusterTree) = typeof(t)

isleaf(clt::ClusterTree)  = isempty(clt.children)
isroot(clt::ClusterTree)  = clt.parent == clt
hasdata(clt::ClusterTree) = isdefined(clt,:data)

diameter(node::ClusterTree)                      = node |> container |> diameter
radius(node::ClusterTree)                        = node |> container |> radius
distance(node1::ClusterTree,node2::ClusterTree)  = distance(container(node1), container(node2))

Base.length(node::ClusterTree) = length(index_range(node))

ambient_dimension(clt::ClusterTree) = ambient_dimension(container(clt))

"""
    ClusterTree(elements,splitter;[copy_elements=true])
    ClusterTree{D}(points,splitter;[copy_elements=true])

Construct a `ClusterTree` from the  given `elements` using the splitting
strategy encoded in `splitter`. If `copy_elements` is set to false, the
`elements` argument are directly stored in the `ClusterTree` and are permuted
during the tree construction.
"""
function ClusterTree{D}(elements,splitter;copy_elements=true) where {D}
    copy_elements && (elements = deepcopy(elements))
    if splitter isa DyadicSplitter
        # make a cube for bounding box for quad/oct trees
        bbox         = HyperRectangle(elements,true)
    else
        bbox         = HyperRectangle(elements)
    end
    n            = length(elements)
    irange     = 1:n
    loc2glob     = collect(irange)
    children     = nothing
    parent       = nothing
    #build the root, then recurse
    root         = ClusterTree{D}(elements,bbox,irange,loc2glob,children,parent)
    _build_cluster_tree!(root,splitter)
    return root
end
ClusterTree(elements,splitter;copy_elements=true) = ClusterTree{Nothing}(elements,splitter;copy_elements)

function _build_cluster_tree!(current_node,splitter)
    if should_split(current_node,splitter)
        children          = split!(current_node,splitter)
        current_node.children = children
        for child in children
            _build_cluster_tree!(child,splitter)
        end
    end
    return current_node
end

"""
    _binary_split!(node::ClusterTree,dir,pos)
    _binary_split!(f,node::ClusterTree)

Split a `ClusterTree` into two, sorting all elements in the process.

Passing a `dir` and `pos` arguments splits the `bounding_box` box of `node`
along direction `dir` at position `pos`, then sorts all points into the
resulting  left/right nodes.

If passed a predicate `f`, each point is sorted
according to whether `f(x)` returns `true` (point sorted on
the left node) or `false` (point sorted on the right node). At the end a minimal
`HyperRectangle` containing all left/right points is created.
"""
function _binary_split!(f::Function,cluster::ClusterTree{T,S,D},
                        buff= Vector{Int}(undef,length(cluster))) where {T,S,D}
    rec        = container(cluster)
    els        = root_elements(cluster)
    l2g        = loc2glob(cluster)
    irange     = index_range(cluster)
    npts_left  = 0
    npts_right = 0
    xl_left = xl_right = high_corner(rec)
    xu_left = xu_right = low_corner(rec)
    n          = length(irange)
    #sort the points into left and right rectangle
    for i in irange
        pt = els[i] |> coords
        if f(pt)
            xl_left = min.(xl_left,pt)
            xu_left = max.(xu_left,pt)
            npts_left += 1
            buff[npts_left] = i
        else
            xl_right = min.(xl_right,pt)
            xu_right = max.(xu_right,pt)
            buff[n-npts_right] = i
            npts_right += 1
        end
    end
    @assert npts_left + npts_right == n "elements lost during split"
    l2g[irange] = l2g[buff]
    els[irange] = els[buff] # reorders the global index set
    # new ranges for cluster
    left_indices      = irange.start:(irange.start)+npts_left-1
    right_indices     = (irange.start+npts_left):irange.stop
    # bounding boxes
    left_rec   = S(xl_left,xu_left)
    right_rec  = S(xl_right,xu_right)
    # create children
    clt1 = ClusterTree{D}(els,left_rec,left_indices,l2g,nothing,cluster)
    clt2 = ClusterTree{D}(els,right_rec,right_indices,l2g,nothing,cluster)
    return clt1, clt2
end

function _binary_split!(cluster::ClusterTree{T,S,D},dir::Int,pos::Number,
                        buff= Vector{Int}(undef,length(cluster))) where {T,S,D}
    rec        = container(cluster)
    els        = root_elements(cluster)
    l2g        = loc2glob(cluster)
    irange     = index_range(cluster)
    npts_left  = 0
    npts_right = 0
    n          = length(irange)
    left_rec, right_rec = split(rec,dir,pos)
    #sort the points into left and right rectangle
    for i in irange
        pt = els[i] |> coords
        if pt in left_rec
            npts_left += 1
            buff[npts_left] = i
        else # pt in right_rec
            buff[n-npts_right] = i
            npts_right += 1
        end
    end
    @assert npts_left + npts_right == n "elements lost during split"
    l2g[irange] = l2g[buff]
    els[irange] = els[buff] # reorders the global index set
    # new ranges for cluster
    left_indices      = irange.start:(irange.start)+npts_left-1
    right_indices     = (irange.start+npts_left):irange.stop
    # create children
    clt1 = ClusterTree{D}(els,left_rec,left_indices,l2g,nothing,cluster)
    clt2 = ClusterTree{D}(els,right_rec,right_indices,l2g,nothing,cluster)
    return clt1, clt2
end

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
function should_split(clt,splitter)
    abstract_method(splitter)
end

"""
    split!(clt::ClusterTree,splitter::AbstractSplitter)

Divide `clt` using the strategy implemented by `splitter`.
"""
function split!(clt,splitter)
    abstract_method(splitter)
end

"""
    struct DyadicSplitter

Used to split an `N` dimensional `ClusterTree` into `2^N` children until at most
`nmax` points are contained in node.
"""
Base.@kwdef struct DyadicSplitter <: AbstractSplitter
    nmax::Int=typemax(Int)
end

function should_split(node::ClusterTree,splitter::DyadicSplitter)
    length(node) > splitter.nmax
end

function split!(cluster::ClusterTree,::DyadicSplitter)
    d        = ambient_dimension(cluster)
    clusters = [cluster]
    rec = cluster.container
    rec_center = center(rec)
    for i in 1:d
        pos = rec_center[i]
        nel = length(clusters) #2^(i-1)
        for _ in 1:nel
            clt = popfirst!(clusters)
            append!(clusters,_binary_split!(clt,i,pos))
        end
    end
    return clusters
end

"""
    struct GeometricSplitter <: AbstractSplitter

Used to split a `ClusterTree` in half along the largest axis.
"""
@Base.kwdef struct GeometricSplitter <: AbstractSplitter
    nmax::Int=50
end

should_split(node::ClusterTree,splitter::GeometricSplitter) = length(node) > splitter.nmax

function split!(cluster::ClusterTree,splitter::GeometricSplitter)
    rec          = cluster.container
    wmax, imax   = findmax(high_corner(rec) - low_corner(rec))
    left_node, right_node = _binary_split!(cluster, imax, low_corner(rec)[imax]+wmax/2)
    return [left_node, right_node]
end

"""
    struct GeometricMinimalSplitter <: AbstractSplitter

Like [`GeometricSplitter`](@ref), but shrinks the children's containters.
"""
@Base.kwdef struct GeometricMinimalSplitter <: AbstractSplitter
    nmax::Int=50
end

should_split(node::ClusterTree,splitter::GeometricMinimalSplitter) = length(node) > splitter.nmax

function split!(cluster::ClusterTree,splitter::GeometricMinimalSplitter)
    rec  = cluster.container
    wmax, imax  = findmax(high_corner(rec) - low_corner(rec))
    mid = low_corner(rec)[imax]+wmax/2
    predicate = (x) -> x[imax] < mid
    left_node,right_node =  _binary_split!(predicate,cluster)
    return [left_node, right_node]
end

"""
    struct PrincipalComponentSplitter <: AbstractSplitter
"""
@Base.kwdef struct PrincipalComponentSplitter <: AbstractSplitter
    nmax::Int=50
end

should_split(node::ClusterTree,splitter::PrincipalComponentSplitter) = length(node) > splitter.nmax

function split!(cluster::ClusterTree,splitter::PrincipalComponentSplitter)
    pts       = cluster._elements
    irange    = cluster.index_range
    xc        = center_of_mass(cluster)
    # compute covariance matrix for principal direction
    cov  = sum(irange) do i
        x = coords(pts[i])
        (x - xc)*transpose(x - xc)
    end
    v = eigvecs(cov)[:,end]
    predicate = (x) -> dot(x-xc,v) < 0
    left_node, right_node = _binary_split!(predicate,cluster)
    return [left_node, right_node]
end

function center_of_mass(clt::ClusterTree)
    pts       = clt._elements
    loc_idxs  = clt.index_range
    # w    = clt.weights
    n    = length(loc_idxs)
    # M    = isempty(w) ? n : sum(i->w[i],glob_idxs)
    # xc   = isempty(w) ? sum(i->pts[i]/M,glob_idxs) : sum(i->w[i]*pts[i]/M,glob_idxs)
    M    = n
    xc   = sum(i->coords(pts[i])/M,loc_idxs)
    return xc
end

"""
    struct CardinalitySplitter <: AbstractSplitter

Used to split a `ClusterTree` along the largest dimension if
`length(tree)>nmax`. The split is performed so the `data` is evenly distributed
amongst all children.
"""
@Base.kwdef struct CardinalitySplitter <: AbstractSplitter
    nmax::Int=50
end

should_split(node::ClusterTree,splitter::CardinalitySplitter) = length(node) > splitter.nmax

function split!(cluster::ClusterTree,splitter::CardinalitySplitter)
    points     = cluster._elements
    irange     = cluster.index_range
    rec        = container(cluster)
    _, imax    = findmax(high_corner(rec) - low_corner(rec))
    med        = median(coords(points[i])[imax] for i in irange) # the median along largest axis `imax`
    predicate = (x) -> x[imax] < med
    left_node, right_node = _binary_split!(predicate,cluster)
    return [left_node, right_node]
end

function Base.show(io::IO,tree::ClusterTree{T,S,D}) where {T,S,D}
    print(io,"ClusterTree with $(length(tree.index_range)) elements of type $T")
end

function Base.summary(clt::ClusterTree)
    @printf "Cluster tree with %i elements" length(clt)
    nodes = collect(AbstractTrees.PreOrderDFS(clt))
    @printf "\n\t number of nodes: %i" length(nodes)
    leaves = collect(AbstractTrees.Leaves(clt))
    @printf "\n\t number of leaves: %i" length(leaves)
    points_per_leaf = map(length,leaves)
    @printf "\n\t min number of elements per leaf: %i" minimum(points_per_leaf)
    @printf "\n\t max number of elements per leaf: %i" maximum(points_per_leaf)
    depth_per_leaf = map(depth,leaves)
    @printf "\n\t min depth of leaves: %i" minimum(depth_per_leaf)
    @printf "\n\t max depth of leaves: %i" maximum(depth_per_leaf)
end
