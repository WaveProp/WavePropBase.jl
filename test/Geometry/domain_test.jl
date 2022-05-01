using Test
import WavePropBase as WPB

# Two points
WPB.clear_entities!()
p1 = WPB.ElementaryEntity(0, 1)
p2 = WPB.ElementaryEntity(0, 2)
@test p1 == p1
@test p1 != p2
points = [p1, p2]
# Three lines
l1 = WPB.ElementaryEntity(1, 1, WPB.ElementaryEntity[p1, p2])
l2 = WPB.ElementaryEntity(1, 2, WPB.ElementaryEntity[p1, p2])
l3 = WPB.ElementaryEntity(1, 3, WPB.ElementaryEntity[p1, p2])
lines = [l1, l2, l3]
# Two surfaces
s1 = WPB.ElementaryEntity(2, 1, WPB.ElementaryEntity[l1, l2])
s2 = WPB.ElementaryEntity(2, 2, WPB.ElementaryEntity[l2, l3])
surfaces = [s1, s2]

@testset "ElementaryEntity" begin
    @test WPB.geometric_dimension(p1) == 0
    @test WPB.tag(p1) == 1
    @test WPB.boundary(l1) == points
    for e in vcat(points, lines, surfaces)
        for b in WPB.boundary(e)
            @test WPB.geometric_dimension(e)-1 == WPB.geometric_dimension(b)
        end
    end
end

# Domains
Ω1 = WPB.Domain(s1)
Ω2 = WPB.Domain(s2)
Ω  = WPB.Domain(surfaces)

@testset "Domain" begin
    @test WPB.entities(Ω) == surfaces
    @test length(Ω) == 2
    @test !isempty(Ω)
    @test Ω2 == setdiff(Ω, Ω1)
    @test Ω  == union(Ω1, Ω2)
    @test s1 in Ω
    @test l1 in Ω
    @test Ω[1] == s1
    @test Ω[end] == s2
    for ω in Ω
        for Γ in WPB.boundary(ω)
            @test WPB.geometric_dimension(ω)-1 == WPB.geometric_dimension(Γ)
        end
    end
    @test Ω1 == intersect(Ω1, Ω)
    @test WPB.Domain(l2) == intersect(Ω1, Ω2)
    @test issubset(Ω1, Ω)
    @test Ω2 == setdiff(Ω, Ω1)
    @test union(WPB.Domain.(lines)...) == WPB.skeleton(Ω)
    @test WPB.Domain(l2) == WPB.internal_boundary(Ω)
    @test union(WPB.Domain(l1), WPB.Domain(l3)) == WPB.external_boundary(Ω)
    @test keys(Ω, 1) == [(1,1), (1,2), (1,3)]
    @test keys(Ω) == [(2,1), (2,2)]
    @test WPB.Domain([s1,s2]) == WPB.Domain([s2,s1])
end
