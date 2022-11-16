using Test
import WavePropBase as WPB
using StaticArrays

@testset "Intersection and composite surfaces" begin
    WPB.clear_entities!()
    l1 = WPB.line((0, 0), (1, 0))
    l2 = WPB.line((1, 0), (1, 1))
    l3 = WPB.line((1, 1), (0, 1))
    l4 = WPB.line((0, 1), (0, 0))
    l5 = WPB.line((1, 1), (1, 2))
    l6 = WPB.line((1, 2), (0, 2))
    l7 = WPB.line((0, 2), (0, 1))
    ent1 = WPB.ElementaryEntity(; boundary=[l1, l2, l3, l4])
    ent2 = WPB.ElementaryEntity(; boundary=[l5, l6, l7, WPB.flip_normal(l3)])
    ent3 = WPB.ElementaryEntity(; boundary=[l1, l2, l5, l6, l7, l4])
    Ω1 = WPB.Domain(ent1)
    Ω2 = WPB.Domain(ent2)
    Ω = WPB.union(Ω1, Ω2)
    # intersect. When entities intersect, return them. When they don't, recurse
    # on the boundaries
    @test intersect(Ω, Ω1) == Ω1
    @test intersect(Ω, Ω2) == Ω2
    Γ12 = intersect(Ω1, Ω2)
    Γ21 = intersect(Ω2, Ω1)
    @test Γ12 == WPB.Domain(l3)
    @test Γ12 == Γ21
    # test that the orientation is different for the entity in Γ12 and Γ21
    @test WPB.normal(Γ12[1], 0.5) == SVector(0, 1)
    @test WPB.normal(Γ21[1], 0.5) == SVector(0, -1)
    # check that skeleton includes the interior boundary l3
    Σ = WPB.skeleton(Ω)
    @test l3 ∈ Σ
    @test l3 ∈ WPB.internal_boundary(Ω)
    @test !(l3 ∈ WPB.external_boundary(Ω))
end

@testset "curves" begin
    f = (x) -> SVector(cos(x[1]), sin(x[1]))
    d = WPB.HyperRectangle(0, 2π)
    ent = WPB.ParametricEntity(f, d)
    s = SVector(rand())
    @test ent(s) == f(s)
    jac = WPB.jacobian(ent, s)
    @test jac[1] ≈ -sin(s[1]) && jac[2] ≈ cos(s[1])
    @test ent(s) ≈ WPB.normal(ent, s)
end

@testset "surfaces" begin
    f = (x) -> SVector(x[1], x[2], sin(x[1]))
    d = WPB.HyperRectangle((0.0, 0.0), (1.0, 1.0))
    ent = WPB.ParametricEntity(f, d)
    s = SVector{2}(rand(2))
    @test ent(s) == f(s)
    @test WPB.jacobian(ent, s) ≈ [1 0; 0 1; cos(s[1]) 0]
end

@testset "Parametric entity tests" begin
    @testset "curves" begin
        f = (x) -> SVector(cos(x[1]), sin(x[1]))
        d = WPB.ReferenceLine()
        ent = WPB.ParametricEntity(f, d)
        s = SVector(rand())
        @test ent(s) == f(s)
        jac = WPB.jacobian(ent, s)
        @test jac[1] ≈ -sin(s[1]) && jac[2] ≈ cos(s[1])
        @test ent(s) ≈ WPB.normal(ent, s)
    end
    @testset "surfaces" begin
        f = (x) -> SVector(x[1], x[2], sin(x[1]))
        d = WPB.ReferenceSquare()
        ent = WPB.ParametricEntity(f, d)
        s = SVector{2}(rand(2))
        @test ent(s) == f(s)
        @test WPB.jacobian(ent, s) ≈ [1 0; 0 1; cos(s[1]) 0]
    end
end
