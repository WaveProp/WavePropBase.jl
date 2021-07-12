using Test
using WavePropBase
using WavePropBase.Geometry
using WavePropBase.Interpolation
using StaticArrays

@testset "Tensor product" begin
    cheb_nodes = (n) -> [cos((2*k+1)/2n * π) for k in 0:n-1]
    # 1d
    nodes = cheb_nodes(10)
    f = x->cos(x)
    vals = map(f,nodes)
    p = TensorLagInterp(nodes,vals)
    xtest = 0.1
    @test f(xtest) ≈ p(xtest)
    # 2d
    nx = 10
    ny = 12
    x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
    y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
    f   = (x) -> cos(x[1]*x[2])
    vals = [f((x,y)) for x in x, y in y]
    p   = TensorLagInterp((x,y),vals)
    xtest = (0.1,0.2)
    @test p(xtest) ≈ f(xtest)
    # 3d
    nx,ny,nz = 10, 12, 8
    x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
    y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
    z   = [0.5+0.5cos((2k-1)*π/2nz) for k in 1:nz] # Chebyshev nodes
    f   = (x) -> cos(x[1]*x[2])*x[3]
    vals = [f((x,y,z)) for x in x, y in y, z in z]
    p   = TensorLagInterp((x,y,z),vals)
    xtest = (0.1,0.2,-0.1)
    @test p(xtest) ≈ f(xtest)
end

# @testset "Cheb tensor product" begin
#     cheb_nodes = (n) -> [-cos(k*π/(n-1)) for k in 0:n-1]
#     # 1d
#     nodes = cheb_nodes(10)
#     f = x->cos(x)
#     vals = map(f,nodes)
#     domain = HyperRectangle(-1,1)
#     p = ChebInterp(vals,domain)
#     @test nodes(p,1) ≈ nodes
#     xtest = 0.1
#     @test f(xtest) ≈ p(xtest)
#     # shift
#     nodes = cheb_nodes(10) .+ 1 .+ π
#     f = x->cos(x)
#     vals = map(f,nodes)
#     domain = HyperRectangle(π,π+2)
#     p = ChebInterp(vals,domain)
#     @test nodes(p,1) ≈ nodes
#     xtest = π+0.1
#     @test f(xtest) ≈ p(xtest)
#     # 2d
# end
