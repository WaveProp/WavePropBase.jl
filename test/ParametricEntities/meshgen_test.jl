using Test
import WavePropBase as WPB
using ForwardDiff

@testset "Circle" begin
    circ = WPB.ParametricEntity(s->(cos(s[1]), sin(s[1])), 0, 2π)
    Γ = WPB.Domain(circ)
    M = WPB.meshgen(Γ, (10,))
    M = WPB.meshgen(Γ; meshsize=0.1)
    # plot(M,Γ)
    @test issetequal(WPB.entities(Γ), [circ])
    @test WPB.geometric_dimension(Γ) == 1
    @test WPB.ambient_dimension(Γ) == 2
end

@testset "Box" begin
    box = WPB.Box() # abstract entity
    Γ = WPB.Domain(WPB.boundary(box))
    M = WPB.meshgen(Γ, (10, 10))
    M = WPB.meshgen(Γ; meshsize=0.1)
    # plot(M,Γ)
    @test WPB.entities(Γ) == WPB.boundary(box)
    @test WPB.geometric_dimension(box) == 3
end

@testset "Ball" begin
    ball = WPB.Ball() # abstract entity
    Γ = WPB.Domain(WPB.boundary(ball))
    M = WPB.meshgen(Γ, (1, 1))
    M = WPB.meshgen(Γ; meshsize=0.1)
    # plot(M,Γ)
    @test WPB.entities(Γ) == WPB.boundary(ball)
    @test WPB.geometric_dimension(ball) == 3
end
