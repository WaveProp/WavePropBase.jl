using Test
using WavePropBase
using WavePropBase.Geometry
using WavePropBase.Interpolation
using StaticArrays

@testset "Generic nodes" begin
    cheb_nodes = (n) -> [cos((2*k+1)/2n * π) for k in 0:n-1]
    # 1d
    nodes = cheb_nodes(10)
    f = x->cos(x)
    vals = map(f,nodes)
    p = TensorLagInterp(vals,nodes)
    xtest = 0.1
    @test f(xtest) ≈ p(xtest)
    # 2d
    nx = 10
    ny = 12
    x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
    y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
    f   = (x) -> cos(x[1]*x[2])
    vals = [f((x,y)) for x in x, y in y]
    p   = TensorLagInterp(vals,(x,y))
    xtest = (0.1,0.2)
    @test p(xtest) ≈ f(xtest)
    # 3d
    nx,ny,nz = 10, 12, 8
    x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
    y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
    z   = [0.5+0.5cos((2k-1)*π/2nz) for k in 1:nz] # Chebyshev nodes
    f   = (x) -> cos(x[1]*x[2])*x[3]
    vals = [f((x,y,z)) for x in x, y in y, z in z]
    p   = TensorLagInterp(vals,(x,y,z))
    xtest = (0.1,0.2,-0.1)
    @test p(xtest) ≈ f(xtest)
end

@testset "Cheb nodes" begin
    # 1d
    n = 10
    vals = zeros(n)
    domain = HyperRectangle(-1,1)
    p      = TensorLagInterp(vals,domain,chebnodes,chebweights)
    f      = x->cos(x[1])
    for I in CartesianIndices(p)
        vals[I] = f(interpolation_nodes(p,I))
    end
    xtest = 0.1
    @test f(xtest) ≈ p(xtest)
    # 2d
    nx = 10
    ny = 12
    vals = zeros(nx,ny)
    domain = HyperRectangle((-2,-1),(0,1))
    p   = TensorLagInterp(vals,domain,chebnodes,chebweights)
    f   = (x) -> cos(x[1]*x[2])
    for I in CartesianIndices(p)
        vals[I] = f(interpolation_nodes(p,I))
    end
    xtest = (-1.,0.2)
    @test p(xtest) ≈ f(xtest)
    # 3d
    nx,ny,nz = 10, 10, 10
    vals = zeros(nx,ny,nz)
    domain = HyperRectangle((-2,-1,2),(0,1,4))
    p   = TensorLagInterp(vals,domain,chebnodes,chebweights)
    f   = (x) -> cos(x[1]*x[2])*sin(x[3])
    for I in CartesianIndices(p)
        vals[I] = f(interpolation_nodes(p,I))
    end
    xtest = (-0.1,0.2,3.1)
    @test p(xtest) ≈ f(xtest)
end
