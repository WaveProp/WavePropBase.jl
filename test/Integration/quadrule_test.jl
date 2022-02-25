using Test
using LinearAlgebra
using WavePropBase
using WavePropBase.Geometry
using WavePropBase.Integration

@testset "Trapezoidal quadrature" begin
    q       = Trapezoidal{10}()
    D       = domain(q)
    @test D == ReferenceLine()
    x,w = q()
    @test sum(w) ≈ 1
    @test all(qnode ∈ D for qnode in x)
    # integrate a periodic function. Should be very accurate.
    @test isapprox(integrate(x->cos(2π*x),q),0,atol=1e-10)
    @test integrate(x->sin(2π*x[1])^2,q) ≈ 0.5
end

@testset "TrapezoidalOpen quadrature" begin
    q = TrapezoidalOpen{10}()
    D = domain(q)
    @test D == ReferenceLine()
    x,w = q()
    @test sum(w) ≈ 1
    @test all(qnode ∈ D for qnode in x)
    # integrate a periodic function. Should be very accurate.
    @test isapprox(integrate(x->cos(2π*x[1]),q),0,atol=1e-10)
    @test integrate(x->sin(2π*x[1])^2,q) ≈ 0.5
end

@testset "Fejer quadrature" begin
    N = 5
    q = Fejer{N}()
    x,w = q()
    D = domain(q)
    @test D == ReferenceLine()
    @test all(qnode ∈ D for qnode in x)
    @test sum(w) ≈ 1
    # integrate all polynomial of degree N-1 exactly
    for n in 1:(N-1)
        @test integrate(x->x[1]^n,q) ≈ 1/(n+1)
    end
end

@testset "Gauss quad on triangle" begin
    d = ReferenceTriangle()
    # exact value for x^a*y^b integrate over reference triangle
    exa = (a,b) -> factorial(a)*factorial(b)/factorial(a+b+2)
    # check all quadrature implemented
    orders = keys(Integration.TRIANGLE_GAUSS_ORDER_TO_NPTS)
    for p in orders
        q     = Gauss(;domain=d,order=p)
        x,w   = q()
        @test domain(q) == d
        @test all(qnode ∈ d for qnode in x)
        for i in 0:p
            for j in 0:p-i
                @test integrate(x->x[1]^i*x[2]^j,q) ≈ exa(i,j)
            end
        end
    end
end

@testset "Gauss quad on tetrahedron" begin
    d = ReferenceTetrahedron()
    for p in (1,2)
        q = Gauss(;domain=d,order=p)
        x,w = q()
        @test sum(w) ≈ 1/6
        @test domain(q) == d
        @test all(qnode ∈ d for qnode in x)
    end
    # FIXME: check that we integrate all monomials up to `order` like in the
    # reference triangle
end

@testset "Tensor product quad on square" begin
    N,M = 10,12
    qx  = Fejer(N)
    qy  = Fejer(M)
    q   = TensorProductQuadrature(qx,qy)
    x,w = q()
    D = domain(q)
    @test D == ReferenceSquare()
    @test all(qnode ∈ D for qnode in x)
    a,b = N-1,M-1 # maximum integration order of monomials
    @test integrate( x -> 1,q) ≈ 1
    f = x -> x[1]^a*x[2]^b
    @test integrate(f,q) ≈ 1/(a+1)*1/(b+1)
end

@testset "Custom quadrature rules" begin
    order = 3
    D = ReferenceTriangle()
    q1 = qrule_for_reference_shape(D,order)
    q2 = CustomQuadratureRule(;domain=D,qnodes=qnodes(q1),qweights=qweights(q1))
    q3 = CustomTriangleQuadratureRule(;qnodes=qnodes(q1),qweights=qweights(q1))
    @test q1() == q2()
    @test q2 == q3
end
