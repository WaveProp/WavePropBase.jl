using Test
import WavePropBase as WPB
using StaticArrays
using LinearAlgebra
using ForwardDiff

@testset "Lagrange elements" begin
    @testset "LagrangeLine" begin
        d = WPB.ReferenceLine()
        f = x -> x[1]^2
        x̂ = WPB.reference_nodes(WPB.LagrangeLine{3})
        vals = f.(x̂)
        p = WPB.LagrangeLine(vals)
        @test p(0) ≈ 0
        @test p(1) ≈ 1
        @test p(0.1) ≈ 0.1^2
        ## line in 3d
        vtx = ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
        l = WPB.LagrangeLine(vtx)
        @test WPB.domain(l) == WPB.ReferenceLine()
        # makes no sense to ask for normal here. Make sure error is thrown.
        @test_throws AssertionError WPB.normal(l, SVector(0.5)) == SVector(-1.0, 1.0) / √2
        ## line in 2d
        a = (0.0, 0.0)
        b = (1.0, 1.0)
        l = WPB.LagrangeLine((a, b))
        @test WPB.domain(l) == WPB.ReferenceLine()
        @test WPB.normal(l, SVector(0.5)) == SVector(1.0, -1.0) / √2
    end
    @testset "LagrangeTriangle" begin
        # triangle in 2d
        vtx = SVector(SVector(0.0, 0.0), SVector(0.0, 1.0), SVector(-1.0, 0))
        t = WPB.LagrangeTriangle(vtx)
        @test WPB.domain_dimension(t) == 2
        @test WPB.range_dimension(t) == 2
        # triangle in 3d
        vtx = SVector(SVector(0.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(-1.0, 0, 0.0))
        t = WPB.LagrangeTriangle(vtx)
        @test WPB.range_dimension(t) == 3
        @test WPB.domain_dimension(t) == 2
        @test WPB.normal(t, SVector(0.1, 0.1)) == SVector(0, 0, 1.0)
        # lagrange nodes in parametric space, following Gmsh convention
        T = SVector{2,Float64}
        triangle3_nodes = [T(0, 0), T(1, 0), T(0, 1)]
        triangle6_nodes = [T(0, 0), T(1, 0), T(0, 1), T(1 / 2, 0), T(1 / 2, 1 / 2),
                           T(0, 1 / 2)]
        triangle10_nodes = [T(0, 0),
                            T(1, 0),
                            T(0, 1),
                            T(1 / 3, 0),
                            T(2 / 3, 0),
                            T(2 / 3, 1 / 3),
                            T(1 / 3, 2 / 3),
                            T(0, 2 / 3),
                            T(0, 1 / 3),
                            T(1 / 3, 1 / 3)]
        triangle15_nodes = [T(0, 0),
                            T(1, 0),
                            T(0, 1),
                            T(0.25, 0),
                            T(0.5, 0),
                            T(0.75, 0),
                            T(0.75, 0.25),
                            T(0.5, 0.5),
                            T(0.25, 0.75),
                            T(0, 0.75),
                            T(0, 0.5),
                            T(0, 0.25),
                            T(0.25, 0.25),
                            T(0.5, 0.25),
                            T(0.25, 0.5)]
        triangle_list = [triangle3_nodes, triangle6_nodes, triangle10_nodes,
                         triangle15_nodes]
        p = T(2 / 5, 3 / 7)   # a point in the Reference Triangle
        @test p ∈ WPB.ReferenceTriangle()
        for triangle_nodes in triangle_list
            Np = length(triangle_nodes)
            # generate points in 3D space
            vtx = @SVector rand(SVector{3,Float64}, Np)
            t = WPB.LagrangeTriangle(vtx)
            # test lagrange basis at lagrange nodes
            for (u, v) in zip(triangle_nodes, vtx)
                @test t(u) ≈ v
            end
            # test derivative
            @test WPB.jacobian(t, p) ≈ ForwardDiff.jacobian(t, p)
        end
    end
    @testset "Tetrahedron" begin
        # TODO: add tests
    end
end

# TODO: move code below to the quadrature tests

# @testset "Line quadrature" begin
#     el     = LagrangeLine((1.,1.),(5.,4.))
#     qrule  = GaussLegendre{1}()
#     x,w    = qrule(el)
#     @test sum(w) ≈ 5
# end

# @testset "Triangle quadrature" begin
#     qrule = Gauss{ReferenceTriangle,1}()
#     F     = LagrangeTriangle((0.,0.),(1.,0),(0.,1.))
#     x,w   = qrule(F)
#     @test sum(w) ≈ 1/2
#     ## equilateral triangle
#     F   = LagrangeTriangle((-1.,0),(1.,0),(0.,1.))
#     x,w = qrule(F)
#     @test sum(w) ≈ 1
# end

# @testset "Triangle surface quadrature" begin
#     qrule = Gauss{ReferenceTriangle,1}()
#     F   = LagrangeTriangle((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))
#     x,w = qrule(F)
#     @test sum(w) ≈ 1/2
#     ## equilateral triangle
#     F   = LagrangeTriangle((-1.,0.,1.),(1.,0.,1.),(0.,1.,1.))
#     x,w = qrule(F)
#     @test sum(w) ≈ 1
# end

# @testset "Tetrahedron quadrature" begin
#     D     = ReferenceTetrahedron
#     qrule = Gauss{D,1}()
#     F   = LagrangeTetrahedron((0,0,0.),(1.,0,0),(0,1,0),(0.,0.,1.))
#     x,w = qrule(F)
#     @test sum(w) ≈ 1/6
#     # dilate by 2x and translate by 1 along  the tetrahedron
#     F   = LagrangeTetrahedron((1,0,0),(3,0,0),(1,2,0),(1,0,2))
#     x,w = qrule(F)
#     @test sum(w) ≈ 1/6*2^3
# end

# @testset "Curved lines" begin
#     qstd = GaussLegendre(20)
#     el   = LagrangeLine((0,0),(1,1),(1/2,1/4))
#     # simple test on smooth integrand
#     f     = x->cos(x[1])
#     Ie     = integrate(f,el)
#     for shand in [IMT(), Kress()]
#         q     = SingularQuadratureRule(qstd,shand)
#         Ia    = integrate(f,q,el)
#         @debug Ia - Ie
#         @test isapprox(Ia,Ie,rtol=1e-6)
#     end
#     # non-smooth integrand
#     f     = x -> log(abs(x[1])*cos(x[1]))
#     Ie     = integrate(f,el)
#     for shand in [IMT(), Kress()]
#         q     = SingularQuadratureRule(qstd,shand)
#         Ia    = integrate(f,q,el)
#         @debug Ia - Ie
#         @test isapprox(Ia,Ie,rtol=1e-5)
#         # check that the `naive` integration woudl have failed the test
#         Istd    = integrate(f,qstd,el)
#         @debug Istd - Ie
#         @test !isapprox(Istd,Ie,rtol=1e-5)
#     end
# end

# @testset "Curved line with singularity inside" begin
#     qstd = GaussLegendre(20)
#     el   = LagrangeLine((0,0),(1,1),(1/2,1/4))
#     vs   = 1/3
#     xs   = el(vs)
#     f    = x -> x == xs ? 0.0 : log(norm(x-xs)*cos(x[1]))
#     Ie    = integrate(f,el)
#     for shand in [IMT(), Kress()]
#         q     = SingularQuadratureRule(qstd,shand)
#         x,w   = q(el,vs)
#         Ia    = integrate(f,x,w)
#         @debug Ie-Ia
#         @test isapprox(Ie,Ia;rtol=1e-5)
#         Ia    = integrate(f,q,el,vs)
#         @test isapprox(Ie,Ia;rtol=1e-5)
#         # check that the `naive` integration woudl have failed the test
#         Istd    = integrate(f,qstd,el)
#         @test !isapprox(Istd,Ie,rtol=1e-5)
#     end
# end
