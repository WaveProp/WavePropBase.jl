using Test
using LinearAlgebra
import WavePropBase as WPB

@testset "Trapezoidal quadrature" begin
    q = WPB.Trapezoidal{10}()
    D = WPB.domain(q)
    @test D == WPB.ReferenceLine()
    x, w = q()
    @test sum(w) ≈ 1
    @test all(qnode ∈ D for qnode in x)
    # integrate a periodic function. Should be very accurate.
    @test isapprox(WPB.integrate(x -> cos(2π * x), q), 0, atol=1e-10)
    @test WPB.integrate(x -> sin(2π * x[1])^2, q) ≈ 0.5
end

@testset "TrapezoidalOpen quadrature" begin
    q = WPB.TrapezoidalOpen{10}()
    D = WPB.domain(q)
    @test D == WPB.ReferenceLine()
    x, w = q()
    @test sum(w) ≈ 1
    @test all(qnode ∈ D for qnode in x)
    # integrate a periodic function. Should be very accurate.
    @test isapprox(WPB.integrate(x -> cos(2π * x[1]), q), 0, atol=1e-10)
    @test WPB.integrate(x -> sin(2π * x[1])^2, q) ≈ 0.5
end

@testset "Fejer quadrature" begin
    N = 5
    q = WPB.Fejer{N}()
    x, w = q()
    D = WPB.domain(q)
    @test D == WPB.ReferenceLine()
    @test all(qnode ∈ D for qnode in x)
    @test sum(w) ≈ 1
    # integrate all polynomial of degree N-1 exactly
    for n in 1:(N - 1)
        @test WPB.integrate(x -> x[1]^n, q) ≈ 1 / (n + 1)
    end
end

@testset "Gauss quad on triangle" begin
    d = WPB.ReferenceTriangle()
    # exact value for x^a*y^b integrate over reference triangle
    exa = (a, b) -> factorial(a) * factorial(b) / factorial(a + b + 2)
    # check all quadrature implemented
    orders = keys(WPB.TRIANGLE_GAUSS_ORDER_TO_NPTS)
    for p in orders
        q = WPB.Gauss(; domain=d, order=p)
        x, w = q()
        @test WPB.domain(q) == d
        @test all(qnode ∈ d for qnode in x)
        for i in 0:p
            for j in 0:(p - i)
                @test WPB.integrate(x -> x[1]^i * x[2]^j, q) ≈ exa(i, j)
            end
        end
    end
end

@testset "Gauss quad on tetrahedron" begin
    d = WPB.ReferenceTetrahedron()
    for p in (1, 2)
        q = WPB.Gauss(; domain=d, order=p)
        x, w = q()
        @test sum(w) ≈ 1 / 6
        @test WPB.domain(q) == d
        @test all(qnode ∈ d for qnode in x)
    end
    # FIXME: check that we integrate all monomials up to `order` like in the
    # reference triangle
end

@testset "Tensor product quad on square" begin
    N, M = 10, 12
    qx = WPB.Fejer(N)
    qy = WPB.Fejer(M)
    q = WPB.TensorProductQuadrature(qx, qy)
    x, w = q()
    D = WPB.domain(q)
    @test D == WPB.ReferenceSquare()
    @test all(qnode ∈ D for qnode in x)
    a, b = N - 1, M - 1 # maximum integration order of monomials
    @test WPB.integrate(x -> 1, q) ≈ 1
    f = x -> x[1]^a * x[2]^b
    @test WPB.integrate(f, q) ≈ 1 / (a + 1) * 1 / (b + 1)
end

@testset "Custom quadrature rules" begin
    order = 3
    D = WPB.ReferenceTriangle()
    q1 = WPB.qrule_for_reference_shape(D, order)
    q2 = WPB.CustomQuadratureRule(;
                                  domain=D,
                                  qcoords=WPB.qcoords(q1),
                                  qweights=WPB.qweights(q1))
    q3 = WPB.CustomTriangleQuadratureRule(;
                                          qcoords=WPB.qcoords(q1),
                                          qweights=WPB.qweights(q1))
    @test q1() == q2()
    @test q2 == q3
end
