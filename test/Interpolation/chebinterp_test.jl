using Test
using WavePropBase
using WavePropBase.Geometry
using WavePropBase.Interpolation
using StaticArrays

@testset "Cheb tensor product" begin
    # 1d
    n = 10
    xcheb = [Interpolation.chebnode(i,n) for i in 1:n]
    f = x->cos(x)*exp(x)
    vals = map(f,xcheb)
    domain = HyperRectangle(-1,1)
    p = ChebInterp(vals,domain)
    @test interpolation_nodes(p,1) ≈ xcheb
    xtest = 0.1
    @test f(xtest) ≈ p(xtest)
    # shift
    xcheb =  xcheb .+ (1 + π)
    f = x->cos(x)*exp(x)
    vals = map(f,xcheb)
    domain = HyperRectangle(π,π+2)
    p = ChebInterp(vals,domain)
    @test interpolation_nodes(p,1) ≈ xcheb
    xtest = π+0.1
    @test f(xtest) ≈ p(xtest)
    # 2d
    # domain = HyperRectangle((-1,0),(2,1))
    # xcheb = cheb_nodes(10) .+ 1 .+ π
end
