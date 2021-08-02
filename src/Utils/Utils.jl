"""
    Utils

Module containing various utility functions for `WaveProp`.
"""
module Utils

using StaticArrays
using WavePropBase

# import all methods in WavePropBase.INTERFACE_METHODS
WavePropBase.@import_interface

export
    # Type alises
    SType,
    # methods
    svector,
    matrix_to_blockmatrix,
    blockmatrix_to_matrix,
    diagonalblockmatrix_to_matrix,
    blockvector_to_vector,
    vector_to_blockvector,
    notimplemented,
    abstractmethod,
    assert_extension,
    assert_concrete_type,
    cross_product_matrix,
    fill_zero_diagonal!,
    Point3D,
    fill_zero_diagonal!,
    cart2sph,
    cart2pol,
    sph2cart,
    pol2cart,
    getnodes

"""
    svector(f,n)

Just like [`Base.ntuple`](https://docs.julialang.org/en/v1/base/base/#Base.ntuple), but convert output to an `SVector`.
"""
@inline svector(f,n) = ntuple(f,n) |> SVector

"""
    blockmatrix_to_matrix(A::Matrix{B}) where {B<:SMatrix}

Convert a `Matrix{B}`, where `B<:SMatrix`, to the equivalent `Matrix{T}`, where `T = eltype(B)`
"""
function blockmatrix_to_matrix(A::Matrix{B})  where B <: SMatrix
    T = eltype(B)
    sblock = size(B)
    ss     = size(A).*sblock # matrix size when viewed as matrix over T
    Afull = Matrix{T}(undef,ss)
    for i=1:ss[1], j=1:ss[2]
        bi, ind_i = divrem(i-1,sblock[1]) .+ (1,1)
        bj, ind_j = divrem(j-1,sblock[2]) .+ (1,1)
        Afull[i,j] = A[bi,bj][ind_i,ind_j]
    end
    return Afull
end

"""
    diagonalblockmatrix_to_matrix(A::Matrix{B}) where {B<:SMatrix}

Convert a diagonal block matrix `A::AbstractVector{B}`, where `A` is the list of diagonal blocks
and `B<:SMatrix`, to the equivalent `SparseMatrixCSC{T}`, where `T = eltype(B)`.
"""
function diagonalblockmatrix_to_matrix(A::AbstractVector{B}) where B<:SMatrix
    T = eltype(B)
    sblock = size(B)
    ss = size(A) .* sblock  # matrix size when viewed as matrix over T
    I = Int64[]
    J = Int64[]
    V = T[]
    i_full, j_full = (1, 1)
    for subA in A
        i_tmp = i_full
        for j in 1:sblock[2]
            i_full = i_tmp
            for i in 1:sblock[1]
                push!(I, i_full)
                push!(J, j_full)
                push!(V, subA[i, j])
                i_full += 1
            end
            j_full += 1
        end
    end
    return sparse(I, J, V, ss[1], ss[2])
end

"""
    blockvector_to_vector(A::Vector{B}) where {B<:SVector}

Convert a `Vector{B}`, where `B<:SVector`, to the equivalent `Vector{T}`, where `T = eltype(B)`
"""
function blockvector_to_vector(A::Vector{B})  where B <: SVector
    T = eltype(B)
    reinterpret(T,A) |> collect
end

"""
    matrix_to_blockmatrix(A::Matrix,B)

Convert a `Matrix{T}` to a `Matrix{B}`, where `B<:Type{SMatrix}`. The element
type of `B` must match that of `A`, and the size of `A` must be divisible by the
size of `B` along each dimension.
"""
function matrix_to_blockmatrix(A::Matrix,B::Type{<:SMatrix})
    @assert eltype(A) == eltype(B)
    @assert sum(size(A) .% size(B)) == 0 "block size $(size(B)) not compatible with size of A=$(size(A))"
    sblock = size(B)
    nblock = div.(size(A),sblock)
    Ablock = Matrix{B}(undef,nblock)
    for i in 1:nblock[1]
        istart = (i-1)*sblock[1] + 1
        iend = i*sblock[1]
        for j in 1:nblock[2]
            jstart = (j-1)*sblock[2] + 1
            jend   = j*sblock[2]
            Ablock[i,j] = A[istart:iend,jstart:jend]
        end
    end
    return Ablock
end

"""
    vector_to_blockvector(A::Vector,B)

Convert a `Vector{T}` to a `Vector{B}`, where `B<:Type{SVector}`. The element
type of `B` must match that of `A`, and the size of `A` must be divisible by the
size of `B` along each dimension.
"""
function vector_to_blockvector(A::Vector,B::Type{<:SVector})
    @assert eltype(A) == eltype(B)
    @assert sum(size(A) .% size(B)) == 0 "block size $(size(B)) not compatible with size of A=$(size(A))"
    reinterpret(B,A) |> collect
end

"""
    notimplemented()

Things which should probably be implemented at some point.
"""
function notimplemented()
    error("not (yet) implemented")
end

"""
    abstractmethod

A method of an `abstract type` for which concrete subtypes are expected
to provide an implementation.
"""
function abstractmethod(T::DataType)
    error("this method needs to be implemented by the concrete subtype $T.")
end
abstractmethod(x) = abstractmethod(typeof(x))

"""
    assert_extension(fname,ext,[msg])

Check that `fname` is of extension `ext`. Print the message `msg` as an assertion error otherwise.
"""
function assert_extension(fname::String,ext::String,msg="file extension must be $(ext)")
    r = Regex("$(ext)\$")
    @assert occursin(r,fname) msg
end

function assert_concrete_type(T::DataType)
    isconcretetype(T) || throw(ConcreteInferenceError(T))
end

"""
    enable_debug(mname)

Activate debugging messages.
"""
function enable_debug()
    ENV["JULIA_DEBUG"] = WavePropBase
end

struct ConcreteInferenceError <: Exception
    T
end

Base.showerror(io::IO, e::ConcreteInferenceError) = print(io, "unable to infer concrete type from function signature: T=$(e.T)" )

"""
    print_threads_info()

Prints in console the total number of threads.
"""
function print_threads_info()
    @info "Number of threads: $(Threads.nthreads())"
end

"""
    const Point1D
    const Point1D(x1)

A point in 1D space, stored in a StaticArray.
Point1D = SVector{1, Float64}.
"""
const Point1D = SVector{1, Float64}

"""
    const Point2D
    const Point2D(x1, x2)
    const Point2D(x::NTuple{2, Float64})

A point in 2D space, stored in a StaticArray.
Point2D = SVector{2, Float64}.
"""
const Point2D = SVector{2, Float64}

"""
    const Point3D
    const Point3D(x1, x2, x3)
    const Point3D(x::NTuple{3, Float64})

A point in 3D space, stored in a StaticArray.
Point3D = SVector{3, Float64}.
"""
const Point3D = SVector{3, Float64}

"""
    const ComplexPoint3D
    const ComplexPoint3D(x1, x2, x3)
    const ComplexPoint3D(x::NTuple{3, ComplexF64})

A complex 3D point, stored in a StaticArray.
ComplexPoint3D = SVector{3, ComplexF64}.
"""
const ComplexPoint3D = SVector{3, ComplexF64}

"""
    const ComplexPoint2D
    const ComplexPoint2D(x1, x2)
    const ComplexPoint2D(x::NTuple{2, ComplexF64})

A complex 2D point, stored in a StaticArray.
ComplexPoint2D = SVector{2, ComplexF64}.
"""
const ComplexPoint2D = SVector{2, ComplexF64}

"""
    const SType{T} = Union{T,Type{T}}

Union type of `T` and its data type `Type{T}`. Used to simplify methods defined
on singleton types where both `foo(::T)` and `foo(::Type{T})` are required.
"""
const SType{T} = Union{T,Type{T}}

"""
    sort_by_type(v)

Sort the elements of `v` into vectors `vi` according to their type. Return a
`Dict{DataType,Vector}` mapping each type to a vector of that type.

# Examples
```julia
v = [1,"a",3,"b"]
dict = sort_by_type(v)
```
"""
function sort_by_type(v)
    dict = Dict{DataType,Vector}()
    for el in v
        T = typeof(el)
        vi = get!(dict,T) do
            T[]
        end
        push!(vi,el)
    end
    return dict
end

"""
    cart2sph(x,y,z)

Map cartesian coordinates `x,y,z` to spherical ones `r, θ, φ` representing the
radius, elevation, and azimuthal angle respectively. The convention followed is
that `0 ≤ θ ≤ π` and ` -π < φ ≤ π`.
"""
function cart2sph(x,y,z)
    azimuth   = atan(y,x)
    elevation = atan(sqrt(x^2 + y^2),z)
    r = sqrt(x^2 + y^2 + z^2)
    return r, elevation, azimuth
end

"""
    sph2cart(x,y,z)

Map spherical coordinates `r,θ,φ` representing the radius, elevation, and
azimuthal angle respectively, to cartesian ones `x, y, z` .
"""
function sph2cart(r,θ,φ)
    x = r*cos(φ)*sin(θ)
    y = r*sin(φ)*sin(θ)
    z = r*cos(θ)
    return x,y,z
end

"""
    cart2pol(x,y)

Map cartesian coordinates `x,y` to polar coordinates `r,θ`. The convention
followed is that `-π ≤ θ ≤ π`.
"""
function cart2pol(x,y)
    r = sqrt(x^2 + y^2)
    θ = atan(y,x)
    return r,θ
end

"""
    pol2cart(r,θ)

Map polar coordinates `r,θ` to cartesian coordinates `x,y`.
"""
function pol2cart(r,θ)
    x = r*cos(θ)
    y = r*sin(θ)
    return x,y
end

"""
    cross_product_matrix(v)

Returns the matrix `Aᵥ` associated with the
cross product `v × ϕ` so that `v × ϕ = Aᵥϕ`.
"""
function cross_product_matrix(v)
    return transpose(SMatrix{3,3,Float64,9}(      0, -v[3],  v[2],
                                             v[3],       0, -v[1],
                                            -v[2],    v[1],     0))
end

function fill_zero_diagonal!(A::AbstractMatrix)
    # CHECK: is this function necessary?
    # it seems that, by default, the operators'
    # diagonals are always zero
    n, m = size(A)
    @assert n == m
    for i in 1:n
        A[i, i] = zero(eltype(A))
    end
end

"""
    @interface f [n=1]

Declare that the function `f` is an interface function. The call
`f(args...)` resolves to `M.f(args...)` where `M` is parent module of the
`args[n]` object.

The somewhat contrived example below illustrates how this can be used to have a
generic method defined in module `A` applied to a type defined on module `B`
which is independent of `A` but which implements the interface
function `f`:
```jldoctest
module A
    using WavePropBase.Utils
    Utils.@interface foo
    # a method which works on any type `x` implementing the `foo` function
    do_work(x) = 2*foo(x)
end

module B
    struct Foo end
    foo(x::Foo) = 1
end

using .A
using .B
foo = B.Foo()
A.do_work(foo)

# output

2
```

Note that if in the example above module `A` implements a generic version of
`foo`, the call `A.do_work(foo)` would use that method instead based on the
dispatch rules.
"""
macro interface(f,n=1)
    ex = quote
        @generated function $f(args...)
            M = parentmodule(args[$n])
            hasmethod(M.$f,args) || error("function `$(M.$f)` must be implemented in module `$M`")
            return quote
                $(M.$f)(args...)
            end
        end
    end
    return esc(ex)
end

"""
    depth(node,acc=0)

Recursive function to compute the depth of `node` in a a tree-like structure.
Require the method `getparent(node)` to be implemented. Overload this function
if your structure has a more efficient way to compute `depth` (e.g. if it stores
it in a field).
"""
function depth(node,acc=0)
    if isroot(node)
        return acc
    else
        depth(node.parent,acc+1)
    end
end

"""
    getnodes(filter,tree,[isterminal=true])

Return all the nodes of `tree` satisfying `filter(node)::Bool`. If
`isterminal`, do not consider childrens of a node for which
`filter(block)==true`.
"""
function getnodes(f,tree,isterminal=true)
    blocks = Vector{typeof(tree)}()
    getnodes!(f,blocks,tree,isterminal)
end

"""
    getnodes!(filter,nodes,tree,[isterminal=true])

Like [`getnodes`](@ref), but append valid nodes to `nodes`.
"""
function getnodes!(f,blocks,tree,isterminal)
    if f(tree)
        push!(blocks,tree)
        # terminate the search along this path if terminal=true
        isterminal || map(x->getnodes!(f,blocks,x,isterminal),getchildren(tree))
    else
        # continue on on children
        map(x->getnodes!(f,blocks,x,isterminal),tree.children)
    end
    return blocks
end

end # module
