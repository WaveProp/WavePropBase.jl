using Test
using WavePropBase
using StaticArrays

function cheb_nodes(n)
    [cos((2*k+1)/2n * π) for k in 0:n-1]
end

@testset "Tensor product" begin
    # 1d
    nodes = cheb_nodes(10)
    f = x->cos(x)
    vals = map(f,nodes)
    p = LagrangeInterp(nodes,vals)
    xtest = 0.1
    @test f(xtest) ≈ p(xtest)
    # 2d
    nx = 10
    ny = 12
    x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
    y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
    f   = (x) -> cos(x[1]*x[2])
    vals = [f((x,y)) for x in x, y in y]
    p   = LagrangeInterp((x,y),vals)
    xtest = (0.1,0.2)
    @test p(xtest) ≈ f(xtest)
    # 3d
    nx,ny,nz = 10, 12, 8
    x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
    y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
    z   = [0.5+0.5cos((2k-1)*π/2nz) for k in 1:nz] # Chebyshev nodes
    f   = (x) -> cos(x[1]*x[2])*x[3]
    vals = [f((x,y,z)) for x in x, y in y, z in z]
    p   = LagrangeInterp((x,y,z),vals)
    xtest = (0.1,0.2,-0.1)
    @test p(xtest) ≈ f(xtest)
end
