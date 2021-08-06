"""
    mutable struct ClusterTree{N,T,D}

Tree structure used to sort points in `N` dimensions. A `data` field of type `D`
can be associated with each node (it defaults to `D=Nothing`).

# Fields:
- `points::Vector{SVector{N,T}}` : vector containing all points present in the
  `root` of the tree. The `points` are permuted during the construction of the
  tree so as to make them contiguous in each node.
- `loc_idxs::UnitRange{Int}` : indices of points contained in the current node.
- `bounding_box::HyperRectangle{N,T}` : axis-aligned bounding box containing all
  points in the node.
- `loc2glob::Vector{Int}` : permutation from global local indexing system to the
  original (global) indexing system used before the construction of the tree.
- `glob2loc::Vector{Int}` : inverse of `loc2glob` permutation.
- `children::Vector{ClusterTree{N,T,D}}`
- `parent::ClusterTree{N,T,D}`
- `data::D` : generic data field parametrically typed on `D`.
"""
mutable struct ClusterTree{N,T,D}
    points::Vector{SVector{N,T}}
    loc_idxs::UnitRange{Int}
    bounding_box::HyperRectangle{N,T}
    loc2glob::Vector{Int}
    glob2loc::Vector{Int}
    children::Vector{ClusterTree{N,T,D}}
    parent::ClusterTree{N,T,D}
    data::D
    # inner constructors handling missing fields
    function ClusterTree{D}(points::Vector{SVector{N,T}},loc_idxs,bounding_box,loc2glob,glob2loc,children,parent) where {N,T,D}
        clt = new{N,T,D}(points,loc_idxs,bounding_box,loc2glob,glob2loc)
        clt.children = isnothing(children) ? Vector{typeof(clt)}() : children
        clt.parent   = isnothing(parent)   ? clt : parent
        clt.data     = D()
        return clt
    end
    # if some data is passed, fill the data field
    function ClusterTree(points::Vector{SVector{N,T}},loc_idxs,bounding_box,loc2glob,glob2loc,children,parent,data::D) where {N,T,D}
        clt = new{N,T,D}(points,loc_idxs,bounding_box,loc2glob,glob2loc)
        clt.children = isnothing(children) ? Vector{typeof(clt)}() : children
        clt.parent   = isnothing(parent)   ? clt : parent
        clt.data = data
        return clt
    end
end
ClusterTree(args...;kwargs...) = ClusterTree{Nothing}(args...;kwargs...)

# interface to AbstractTrees. No children is determined by an empty tuple for AbstractTrees
AbstractTrees.children(t::ClusterTree) = isleaf(t) ? () : t.children
AbstractTrees.nodetype(t::ClusterTree) = typeof(t)

isleaf(clt::ClusterTree) = isempty(clt.children)
isroot(clt::ClusterTree) = clt.parent == clt
hasdata(clt::ClusterTree) = isdefined(clt,:data)

children(clt::ClusterTree) = clt.children

diameter(node::ClusterTree)                      = diameter(node.bounding_box)
radius(node::ClusterTree)                        = diameter(node)/2
distance(node1::ClusterTree,node2::ClusterTree)  = distance(node1.bounding_box, node2.bounding_box)

ambient_dimension(clt::ClusterTree{N}) where {N} = N

Base.length(node::ClusterTree) = length(node.loc_idxs)
Base.range(node::ClusterTree)  = node.loc_idxs

"""
    ClusterTree(points,splitter;[copy_points=true])
    ClusterTree{D}(points,splitter;[copy_points=true])

Construct a `ClusterTree` from the  given `points` using the splitting strategy
encoded in `splitter`. If `copy_points` is set to false, the `points` argument
is directly stored in the `ClusterTree`, and is therefore permuted into the
local indexing.
"""
function ClusterTree{N,T,D}(points::Vector{SVector{N,T}},splitter;copy_points=true) where {N,T,D}
    copy_points && (points = deepcopy(points))
    bbox         = HyperRectangle(points)
    n            = length(points)
    loc_idxs     = 1:n
    loc2glob     = collect(loc_idxs)
    glob2loc     = copy(loc2glob)
    children     = nothing
    parent       = nothing
    #build the root, then recurse
    root         = ClusterTree{D}(points,loc_idxs,bbox,loc2glob,glob2loc,children,parent)
    _build_cluster_tree!(root,splitter)
    root.glob2loc = invperm(root.loc2glob)
    permute!(points,glob2loc)
    return root
end
ClusterTree{D}(points::Vector{SVector{N,T}},spl) where {N,T,D} =  ClusterTree{N,T,D}(points,spl)

function _build_cluster_tree!(current_node,splitter)
    if should_split(current_node,splitter)
        children          = split!(current_node,splitter)
        current_node.children = children
        for child in children
            child.parent = current_node
            _build_cluster_tree!(child,splitter)
        end
    end
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
`nmax` points are contained in node *or* the depth `dmax` is reached.
"""
Base.@kwdef struct DyadicSplitter <: AbstractSplitter
    nmax::Int=typemax(Int)
    dmax::Int=-1
end

function should_split(node::ClusterTree,splitter::DyadicSplitter)
    length(node) > splitter.nmax || depth(node) < splitter.dmax
end

function split!(cluster::ClusterTree,splitter::DyadicSplitter)
    d        = ambient_dimension(cluster)
    clusters = [cluster]
    for i in 1:d
        rec  = cluster.bounding_box
        pos = (rec.high_corner[i] + rec.low_corner[i])/2
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
    rec          = cluster.bounding_box
    wmax, imax   = findmax(rec.high_corner - rec.low_corner)
    left_node, right_node = _binary_split!(cluster, imax, rec.low_corner[imax]+wmax/2)
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
    rec  = cluster.bounding_box
    wmax, imax  = findmax(rec.high_corner - rec.low_corner)
    mid = rec.low_corner[imax]+wmax/2
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
    pts       = cluster.points
    loc_idxs  = cluster.loc_idxs
    glob_idxs = view(cluster.loc2glob,loc_idxs)
    xc   = center_of_mass(cluster)
    cov  = sum(glob_idxs) do i
        (pts[i] - xc)*transpose(pts[i] - xc)
    end
    v = eigvecs(cov)[:,end]
    predicate = (x) -> dot(x-xc,v) < 0
    left_node, right_node = _binary_split!(predicate,cluster)
    return [left_node, right_node]
end

function center_of_mass(clt::ClusterTree)
    pts       = clt.points
    loc_idxs  = clt.loc_idxs
    glob_idxs = view(clt.loc2glob,loc_idxs)
    # w    = clt.weights
    n    = length(loc_idxs)
    # M    = isempty(w) ? n : sum(i->w[i],glob_idxs)
    # xc   = isempty(w) ? sum(i->pts[i]/M,glob_idxs) : sum(i->w[i]*pts[i]/M,glob_idxs)
    M    = n
    xc   = sum(i->pts[i]/M,glob_idxs)
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
    points     = cluster.points
    loc_idxs   = cluster.loc_idxs
    glob_idxs  = view(cluster.loc2glob,loc_idxs)
    rec        = cluster.bounding_box
    _, imax  = findmax(rec.high_corner - rec.low_corner)
    med         = median(points[i][imax] for i in glob_idxs) # the median along largest axis `imax`
    predicate = (x) -> x[imax] < med
    left_node, right_node = _binary_split!(predicate,cluster)
    return [left_node, right_node]
end

"""
    _binary_split!(node::ClusterTree,dir,pos)
    _binary_split!(f,node::ClusterTree)

Split a `ClusterTree` into two, sorting all points and data in the process.

Passing a `dir` and `pos` arguments splits the `bounding_box` box of `node`
along direction `dir` at position `pos`, then sorts all points into the
resulting  left/right nodes.

If passed a predicate `f`, each point is sorted
according to whether `f(x)` returns `true` (point sorted on
the left node) or `false` (point sorted on the right node). At the end a minimal
`HyperRectangle` containing all left/right points is created.
"""
function _binary_split!(cluster::ClusterTree{N,T,D},dir::Int,pos::Number) where {N,T,D}
    points        = cluster.points
    # weights       = cluster.weights
    loc_idxs      = cluster.loc_idxs
    glob_idxs     = view(cluster.loc2glob,loc_idxs)
    glob_idxs_new = view(cluster.glob2loc,loc_idxs)
    npts_left  = 0
    npts_right = 0
    rec = cluster.bounding_box
    left_rec, right_rec = split(rec,dir,pos)
    n                   = length(loc_idxs)
    #sort the points into left and right rectangle
    for i in glob_idxs
        pt = points[i]
        if pt in left_rec
            npts_left += 1
            glob_idxs_new[npts_left]    = i
        else  pt # pt in right_rec
            glob_idxs_new[n-npts_right] = i
            npts_right += 1
        end
    end
    @assert npts_left + npts_right == n "points lost during split"
    # update loc2glob map
    copy!(glob_idxs,glob_idxs_new)
    # new ranges for cluster
    left_indices      = loc_idxs.start:(loc_idxs.start)+npts_left-1
    right_indices     = (loc_idxs.start+npts_left):loc_idxs.stop
    # create children
    clt1 = ClusterTree{D}(points,left_indices,  left_rec,  cluster.loc2glob, cluster.glob2loc,nothing, cluster)
    clt2 = ClusterTree{D}(points,right_indices, right_rec, cluster.loc2glob, cluster.glob2loc,nothing, cluster)
    return clt1, clt2
end

function _binary_split!(f::Function,cluster::ClusterTree{N,T,D}) where {N,T,D}
    points        = cluster.points
    # weights       = cluster.weights
    loc_idxs      = cluster.loc_idxs
    glob_idxs     = view(cluster.loc2glob,loc_idxs)
    glob_idxs_new = view(cluster.glob2loc,loc_idxs)
    npts_left  = 0
    npts_right = 0
    xl_left = xl_right = svector(i->typemax(T),N)
    xu_left = xu_right = svector(i->typemin(T),N)
    n          = length(loc_idxs)
    #sort the points into left and right rectangle
    for i in glob_idxs
        pt = points[i]
        if f(pt)
            xl_left = min.(xl_left,pt)
            xu_left = max.(xu_left,pt)
            npts_left += 1
            glob_idxs_new[npts_left]    = i
        else
            xl_right = min.(xl_right,pt)
            xu_right = max.(xu_right,pt)
            glob_idxs_new[n-npts_right] = i
            npts_right += 1
        end
    end
    @assert npts_left + npts_right == n "points lost during split"
    # update loc2glob map
    copy!(glob_idxs,glob_idxs_new)
    # new ranges for cluster
    left_indices      = loc_idxs.start:(loc_idxs.start)+npts_left-1
    right_indices     = (loc_idxs.start+npts_left):loc_idxs.stop
    # compute bounding boxes
    left_rec   = HyperRectangle(xl_left,xu_left)
    right_rec  = HyperRectangle(xl_right,xu_right)
    # create children
    clt1 = ClusterTree{D}(points,left_indices,left_rec,cluster.loc2glob, cluster.glob2loc,nothing,cluster)
    clt2 = ClusterTree{D}(points,right_indices,right_rec,cluster.loc2glob, cluster.glob2loc,nothing,cluster)
    return clt1, clt2
end

function Base.show(io::IO,tree::ClusterTree{N,T,D}) where {N,T,D}
    print(io,"ClusterTree{$N,$T,$D} with $(length(tree)) points")
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
