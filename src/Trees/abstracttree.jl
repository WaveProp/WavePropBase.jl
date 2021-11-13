"""
    abstract type AbstracTree end

Supertype for tree-like objects.
"""
abstract type AbstractTree end

"""
    children(t::AbstractTree)

Iterable collection of the node's children.
"""
function children(t::AbstractTree)
    abstractmethod(t)
end

"""
    parent(t::AbstractTree)

The node's parent. If `t` is a root, then `parent(t)==t`.
"""
function Base.parent(t::AbstractTree)
    abstractmethod(t)
end

"""
    isroot(t)

Return `true` if `t` is its own parent.
"""
function isroot(t)
    parent(t) == t
end

"""
    isleaf(t)

Return `true` if `t` has no children.
"""
function isleaf(t::AbstractTree)
    isempty(children(t))
end

"""
    depth(tree::AbstractTree,acc=0)

Recursive function to compute the depth of `node` in a a tree-like structure.

Overload this function if your structure has a more efficient way to compute
`depth` (e.g. if it stores it in a field).
"""
function depth(tree,acc=0)
    if isroot(tree)
        return acc
    else
        depth(parent(tree),acc+1)
    end
end

"""
    filter_tree(f,tree,isterminal=true)

Return a vector containing all the nodes of `tree` such that
`filter(node)==true`.  The argument `isterminal` can be used to control whether
to continue the search on `children` of nodes for which `f(node)==true`.
"""
function filter_tree(f,tree,isterminal=true)
    nodes = Vector{typeof(tree)}()
    filter_tree!(f,nodes,tree,isterminal)
end

"""
    filter_tree!(filter,nodes,tree,[isterminal=true])

Like [`filter_tree`](@ref), but appends results to `nodes`.
"""
function filter_tree!(f,nodes,tree,isterminal=true)
    if f(tree)
        push!(nodes,tree)
        # terminate the search along this path if terminal=true
        isterminal || map(x->filter_tree!(f,nodes,x,isterminal),children(tree))
    else
        # continue on on children
        map(x->filter_tree!(f,nodes,x,isterminal),children(tree))
    end
    return nodes
end

# interface to AbstractTrees. No children is determined by an empty tuple for
# AbstractTrees.
AbstractTrees.children(t::AbstractTree) = isleaf(t) ? () : t.children
AbstractTrees.nodetype(t::AbstractTree) = typeof(t)
