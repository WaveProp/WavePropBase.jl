```@meta
CurrentModule = WavePropBase
```

# Helmholtz scattering

!!! note "Important points covered in this tutorial"
    - Creating a geometry using [`ParametricEntity`](@ref)
    - Assembling integral operators and integral potential
    - Setting up a sound-soft and transmission problem
    - Solving a boundary integral equation

In this tutorial we will show how to solve an acoustic scattering problem in the
context of Helmholtz equation. We will begin with a simple example involving a
*smooth* sound-soft obstacle, and gradually make the problem more complex
considering a transmission problem. Along the way we will introduce the necessary
techniques used to handle some difficulties encountered. Assuming you
have already installed `WavePropBase`, let us begin by importing the necessary
packages for this tutorial:

```@example helmholtz_scattering_2d
# Load required packages
import WavePropBase as WPB
using LinearAlgebra
using StaticArrays
using Plots
```

We are now ready to begin solving some PDEs!

## Sound-soft scattering

The first example is that of a sound-soft acoustic scattering problem.
Mathematically, this means we will consider an exterior Helmholtz equation
(time-harmonic acoustics) with a Dirichlet boundary condition. More precisely,
let ``\Omega \subset \mathbb{R}^2`` be a bounded domain, and denote by ``\Gamma
= \partial \Omega`` its boundary. Then we wish to solve

```math
    \Delta u + k^2 u = 0 \quad \text{on} \quad \mathbb{R}^2 \setminus \bar{\Omega},
```

subject to Dirichlet boundary conditions on ``\Gamma``:

```math
    u(\boldsymbol{x}) = g(\boldsymbol{x}) \quad \text{for} \quad \boldsymbol{x} \in \Gamma.
```

Here ``g`` is a give data, and ``k`` is a constant.

For concreteness, we will take ``\Gamma`` to be a kite-like object, and focus on
the *plane-wave scattering* problem. This means we will seek a solution ``u`` of
the form ``u = u_s + u_i``, where ``u_i`` is a known incident field, and ``u_s``
is the scattered field we wish to compute. Using the theory of boundary integral
equations, we can express ``u_s`` as

```math
    u_s(\boldsymbol{r}) = \mathcal{D}[\sigma](\boldsymbol{r}) - i k \mathcal{S}[\sigma](\boldsymbol{r}),
```

where ``\mathcal{S}`` is the so-called single layer potential, ``\mathcal{D}``
is the double-layer potential, and ``\sigma : \Gamma \to \mathbb{C}`` is a
surface density. This is an indirect formulation (because ``\sigma`` is an
*auxiliary* density, not necessarily physical) commonly referred to as a
*combined field formulation*. Taking the limit ``\mathbb{R}^2 \setminus \bar
\Omega \ni x \to \Gamma``, it can be shown that the following equation holds on
``\Gamma``: 

```math
    \left( \frac{\mathrm{I}}{2} + \mathrm{D} - i k \mathrm{S} \right)[\sigma] = g,
```

where $\mathrm{I}$ is the identity operator, and $\mathrm{S}$ and $\mathrm{D}$ are the single- and double-layer operators. This is the **combined field integral equation** that we will solve next.

As for the boundary data ``g``, using the sound-soft condition (i.e. ``u=0`` on
the scatterer), it follows that ``u_s = -u_i`` on ``\Gamma``. We are now in a
position to solve the problem! Let us define the **geometry** and the **PDE**
first:

```@example helmholtz_scattering_2d
kite = WPB.ParametricEntity(0,2π) do s
    (cos(s[1]) + 0.65 * cos(2 * s[1]) - 0.65, 1.5 * sin(s[1]))
end
Γ   = WPB.Domain(kite)
k   = 10
pde = WPB.Helmholtz(;dim=2,k)
```

!!! tip "Geometry creation"
    For simple shapes, it is often convenient to define the geometry using a
    `ParametricEntity` as done above. As things get more complex, it is
    recommended to use *gmsh* to handle both the geometrical entities and the
    mesh generation.

Next we need to generate a mesh for ``\Gamma``. This is done using the
[`meshgen`](@ref) function:

```@example helmholtz_scattering_2d
M = WPB.meshgen(Γ;meshsize=0.2)
```

Because we will be employing a Nystrom method, we must append to our mesh
structure information related to the quadrature nodes, normals, etc. This is
done by calling the [`NystromMesh`](@ref) constructor:

```@example helmholtz_scattering_2d
Q = WPB.NystromMesh(M;qorder=3)
```

Let us plot the mesh and quadrature nodes in `Q`

```@example helmholtz_scattering_2d
plot(Q;lw=2,ms=2)
```

!!! note 
    In this simple example, the mesh elements are create by a uniform partition
    of the parameter space.

With the `NystromMesh` constructed, we now can define the integral operators
``\mathrm{S}`` and ``\mathrm{D}`` over ``\Gamma``:

```@example helmholtz_scattering_2d
Sop = WPB.SingleLayerOperator(pde,Q)
Dop = WPB.DoubleLayerOperator(pde,Q)
```

Both `Sop` and `Dop` are of type [`IntegralOperator`](@ref), which is a subtype of
`AbstractMatrix`. They represent a discrete approximation to linear operator
mapping densities defined on the `source_mesh` into densities defined on the
`target_mesh`. There are two well-known difficulties related to the discretization
of these `IntegralOperator`s:
- The kernel of the integral operator is not smooth, and thus specialized
quadrature rules are required to accurately approximate the matrix entries for
which the target and source point lie *close* (relative to some scale) to each
other. 
- The underlying matrix is dense, and thus the storage and computational
cost of the operator is prohibitive for large problems unless acceleration
techniques such as *Fast Multipole Methods* or *Hierarchical Matrices* are
employed.  

In this example, we will simply use a sparse correction method based on the
*density interpolation technique*, and ignore the second issue. More precisely,
because the problem is two-dimensional and simple, we will just assemble a dense
`Matrix` to represent the integral operator.

We first build the *dense* part of the operators:

```@example helmholtz_scattering_2d
S₀ = Matrix(Sop)
D₀ = Matrix(Dop)
```

Next we build the sparse corrections

```@example helmholtz_scattering_2d
δS, δD = WPB.dim_correction(pde,Q,Q,S₀,D₀)
```

We can now add the corrections to the dense part to obtain the final operators:

```@example helmholtz_scattering_2d
S = S₀ + δS
D = D₀ + δD
```

Finally, we are ready to solve the scattering problem:

```@example helmholtz_scattering_2d
# the linear operator
L = 0.5*I + D - im*k*S
# incident wave
θ  = π/4
d⃗  = (cos(θ),sin(θ))
uᵢ = x -> exp(im*k*dot(d⃗,x))
rhs = [-uᵢ(q.coords) for q in Q.qnodes]
σ = L \ rhs
```

We can now reconstruct the solution using the `IntegralPotential` representation:

```@example helmholtz_scattering_2d
𝒮 = WPB.SingleLayerPotential(pde,Q)
𝒟 = WPB.DoubleLayerPotential(pde,Q)
u = x -> 𝒟[σ](x) - im*k*𝒮[σ](x) + uᵢ(x)
```

Note that, albeit not very efficient, the function `u` can be evaluated at any
point. To see what the solution looks like, lets plot it on a grid

```@example helmholtz_scattering_2d
h = 2π/k/10 # 10 pts per wavelength
x = y = -5:h:5
U = [WPB.isinside((x,y),Q) ? NaN + NaN*im : u((x,y)) for y in y, x in x]
fig = heatmap(x,y,real(U),clims=(-2,2))
plot!(fig,M,lc=:black,lw=4)
```







