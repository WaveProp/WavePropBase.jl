using Test
import WavePropBase as WPB
using StaticArrays
using LinearAlgebra
using Random
Random.seed!(1)

@testset "Reference segment" begin
    qstd = WPB.Fejer(20)
    d     = WPB.domain(qstd)
    # simple test on smooth integrand
    for shand in [WPB.IMT(), WPB.Kress(), WPB.KressP()]
        q     = WPB.SingularQuadratureRule(qstd,shand)
        @test isapprox(WPB.integrate(cos,q),sin(1),rtol=1e-2)
    end
    # non-smooth integrand, singularity at 0
    f     = (x) -> log(abs(x))
    Ie    = -1
    for shand in [WPB.IMT(), WPB.Kress(),WPB.KressP()]
        q     = WPB.SingularQuadratureRule(qstd,shand)
        Ia    = WPB.integrate(f,q)
        @test isapprox(Ia,Ie,rtol=1e-4)
        # check that the `naive` integration woudl have failed the test
        Istd    = WPB.integrate(f,qstd)
        @test !isapprox(Istd,Ie,rtol=1e-4)
    end
    # test that KressP() can handle singularity at both endpoints, and Kress() cannot
    f     = (x) -> log(x) + log(1-x)
    Ie    = -2
    q     = WPB.SingularQuadratureRule(qstd,WPB.KressP())
    Ia    = WPB.integrate(f,q)
    @test isapprox(Ia,Ie,rtol=1e-5)
    q     = WPB.SingularQuadratureRule(qstd,WPB.Kress())
    Ia    = WPB.integrate(f,q)
    @test !isapprox(Ia,Ie,rtol=1e-5)
    # non-smooth integrand with singularity inside
    s     = 0.1
    f     = (x) -> log(abs(x-s))
    Ie    = -1.32508297337145
    # FIXME: for the time being, usign IMT with a singular quadrature rule for
    # interior points yields a rounding problem because the change of variables
    # accumulates points very quickly around zero. That is OK when the
    # singularity is at zero, but shifting such nodes to another point s ∈ (0,1]
    # has to be done with care since these can get rounded to `s`, which in turn
    # yields a singular value for the integrand. Note that the same can happen
    # with other singularity handler types, but IMT is worse due to the
    # exponential nature of the change of variables.
    for shand in [WPB.Kress(),WPB.KressP()]
        q     = WPB.SingularQuadratureRule(qstd,shand)
        x,w  = q(s) # singular nodes and weights which depend on s
        Ia    = WPB.integrate(f,x,w)
        @test isapprox(Ia,Ie,rtol=1e-5)
        # check that the `naive` integration woudl have failed the test
        Istd    = WPB.integrate(f,qstd)
        @test !isapprox(Istd,Ie,rtol=1e-5)
    end
end

@testset "Duffy" begin
    q1d   = WPB.Fejer(5)
    qstd  = WPB.TensorProductQuadrature(q1d,q1d)
    duffy = WPB.Duffy()
    qsin  = WPB.SingularQuadratureRule(qstd,duffy)
    x,w  = qsin()
    # regular kernel
    k    = x -> cos(x[1])*sin(x[2])
    Ie   = 1/2*(sin(1)-cos(1))
    Ia   = WPB.integrate(k,qsin)
    @test isapprox(Ie,Ia,rtol=1e-4)
    # singular kernel at right vertex
    s    = SVector(1,0)
    k    = x -> 1/norm(x-s)
    Ia   = WPB.integrate(k,qsin)
    Ie   = acosh(sqrt(2))
    @test isapprox(Ie,Ia,rtol=1e-5)
    # singular kernel at wrong vertex (should not be accurate, see test below)
    s    = SVector(0,1)
    k    = x -> 1/norm(x-s)
    Ie   = acosh(sqrt(2))
    Ia   = WPB.integrate(k,qsin)
    @test !isapprox(Ie,Ia,rtol=1e-5)
end

@testset "2d Kress" begin
    q1d   = WPB.Fejer{10}()
    qstd  = WPB.TensorProductQuadrature(q1d,q1d)
    s1d   = WPB.Kress(order=2)
    sing_handler = WPB.TensorProductSingularityHandler(s1d,s1d)
    qsin = WPB.SingularQuadratureRule(qstd,sing_handler)
    x,w  = qsin()
    # regular kernel
    k = (x) -> cos(x[1])
    Ie   = sin(1)
    Isin = sum(k.(x).*w)
    @test isapprox(Ie,Isin,rtol=1e-6)
    # singular kernel at origin
    k  = (x) -> 1/norm(x)
    Ie = 1/2*log(17 + 12*sqrt(2)) # Wolfram alpha is your friend ;-)
    Isin = sum(k.(x).*w)
    @test isapprox(Ie,Isin,rtol=1e-2)
end
