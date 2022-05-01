using Test
import WavePropBase as WPB
using StaticArrays

@testset "Generic nodes" begin
    cheb_nodes = (n) -> [cos((2*k+1)/2n * π) for k in 0:n-1]
    # 1d
    nodes = cheb_nodes(10)
    f = x->cos(x)
    vals = map(f,nodes)
    p = WPB.TensorLagInterp(vals,nodes)
    xtest = 0.1
    @test f(xtest) ≈ p(xtest)
    # 2d
    nx = 10
    ny = 12
    x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
    y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
    f   = (x) -> cos(x[1]*x[2])
    vals = [f((x,y)) for x in x, y in y]
    p   = WPB.TensorLagInterp(vals,(x,y))
    xtest = (0.1,0.2)
    @test p(xtest) ≈ f(xtest)
    # 3d
    nx,ny,nz = 10, 12, 8
    x   = [0.5+0.5cos((2k-1)*π/2nx) for k in 1:nx] # Chebyshev nodes
    y   = [0.5+0.5cos((2k-1)*π/2ny) for k in 1:ny] # Chebyshev nodes
    z   = [0.5+0.5cos((2k-1)*π/2nz) for k in 1:nz] # Chebyshev nodes
    f   = (x) -> cos(x[1]*x[2])*x[3]
    vals = [f((x,y,z)) for x in x, y in y, z in z]
    p   = WPB.TensorLagInterp(vals,(x,y,z))
    xtest = (0.1,0.2,-0.1)
    @test p(xtest) ≈ f(xtest)
end

@testset "Cheb nodes" begin
    # 1,2,3 dimensions
    for d in 1:3
        p = ntuple(i->10+i,d)
        vals   = zeros(p)
        nodes1d   = map(p->WPB.cheb1nodes(p,-1,1),p)
        weights1d = map(p->WPB.cheb1weights(p),p)
        poly    = WPB.TensorLagInterp(vals,nodes1d,weights1d)
        f     = x->cos(prod(x))
        for I in CartesianIndices(poly)
            vals[I] = f(WPB.interpolation_nodes(poly,I))
        end
        xtest = ntuple(i->0.1,d)
        @test f(xtest) ≈ poly(xtest)
    end
    # 1,2,3 dimensions
    for d in 1:3
        p = ntuple(i->10+i,d)
        vals   = zeros(p)
        nodes1d   = map(p->WPB.cheb2nodes(p,-1.3,1.5),p)
        weights1d = map(p->WPB.cheb2weights(p),p)
        poly    = WPB.TensorLagInterp(vals,nodes1d,weights1d)
        f     = x->cos(prod(x))
        for I in CartesianIndices(poly)
            vals[I] = f(WPB.interpolation_nodes(poly,I))
        end
        xtest = ntuple(i->0.1,d)
        @test f(xtest) ≈ poly(xtest)
        # todo add test when xtest is one of the interpolation points
    end
end
