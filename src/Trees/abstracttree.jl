"""
    abstract type AbstracTree end

Supertype for tree-like objects.
"""
abstract type AbstractTree end

"""
    depth(t::AbstractTree)

Number of predecessors between `t` and the root of the tree.
"""
function depth(t::AbstractTree)
    abstractmethod(t)
end

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
function parent(t::AbstractTree)
    abstractmethod(t)
end

"""
    isroot(t)

Return `true` if `t` is its own parent.
"""
function isroot(t::AbstractTree)
    abstractmethod(t)
end

"""
    isleaf(t)

Return `true` if `t` has no children.
"""
function isleaf(t::AbstractTree)
    abstractmethod(t)
end

"""
    filter(filter,tree::AbstractTree,[isterminal=true])

Return all the nodes of `tree` satisfying `filter(node)==true`. If
`isterminal`, do not recurse on children of nodes for which `filter(node)==true`.
"""
function Base.filter(f,tree::AbstractTree,isterminal=true)
    nodes = Vector{typeof(tree)}()
    filter!(f,nodes,tree,isterminal)
end

"""
    filter!(filter,nodes,tree,[isterminal=true])

Like [`filter`](@ref), but appends results to `nodes`.
"""
function Base.filter!(f,nodes,tree::AbstractTree,isterminal=true)
    if f(tree)
        push!(nodes,tree)
        # terminate the search along this path if terminal=true
        isterminal || map(x->filter!(f,nodes,x,isterminal),getchildren(tree))
    else
        # continue on on children
        map(x->filter!(f,nodes,x,isterminal),tree.children)
    end
    return nodes
end
