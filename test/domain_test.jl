using Test
using WavePropBase

const WPBase = WavePropBase

# Two points
WPBase.clear!()
p1 = WPBase.ElementaryEntity(0, 1, WPBase.ElementaryEntity[])
p2 = WPBase.ElementaryEntity(0, 2, WPBase.ElementaryEntity[])
points = [p1, p2]
# Three lines
l1 = WPBase.ElementaryEntity(1, 1, WPBase.ElementaryEntity[p1, p2])
l2 = WPBase.ElementaryEntity(1, 2, WPBase.ElementaryEntity[p1, p2])
l3 = WPBase.ElementaryEntity(1, 3, WPBase.ElementaryEntity[p1, p2])
lines = [l1, l2, l3]
# Two surfaces
s1 = WPBase.ElementaryEntity(2, 1, WPBase.ElementaryEntity[l1, l2])
s2 = WPBase.ElementaryEntity(2, 2, WPBase.ElementaryEntity[l2, l3])
surfaces = [s1, s2]

@testset "ElementaryEntity" begin
    @test WPBase.key(p1) == (0,1)
    @test WPBase.geometric_dimension(p1) == 0
    @test WPBase.boundary(l1) == points
    for e in vcat(points, lines, surfaces)
        for b in WPBase.boundary(e)
            @test WPBase.geometric_dimension(e)-1 == WPBase.geometric_dimension(b)
        end
    end
end

# Domains
Ω1 = WPBase.Domain(s1)
Ω2 = WPBase.Domain(s2)
Ω  = WPBase.Domain(surfaces)

@testset "Domain" begin
    @test WPBase.entities(Ω) == surfaces
    @test length(Ω) == 2
    @test !isempty(Ω)
    @test Ω2 == setdiff(Ω, Ω1)
    @test Ω == union(Ω1, Ω2)
    @test s1 in Ω
    @test Ω[1] == s1
    @test Ω[end] == s2
    for ω in Ω
        for Γ in WPBase.boundary(ω)
            @test WPBase.geometric_dimension(ω)-1 == WPBase.geometric_dimension(Γ)
        end
    end
    @test Ω1 == intersect(Ω1, Ω)
    @test WPBase.Domain(l2) == intersect(Ω1, Ω2)
    @test issubset(Ω1, Ω)
    @test Ω2 == WPBase.remove(Ω1, Ω)
    @test union(WPBase.Domain.(lines)...) == WPBase.skeleton(Ω)
    @test WPBase.Domain(l2) == WPBase.internal_boundary(Ω)
    @test union(WPBase.Domain(l1), WPBase.Domain(l3)) == WPBase.external_boundary(Ω)
    @test keys(Ω, 1) == [(1,1), (1,2), (1,3)]
    @test keys(Ω) == [(2,1), (2,2)]
    @test WPBase.Domain([s1,s2]) == WPBase.Domain([s2,s1])
end
