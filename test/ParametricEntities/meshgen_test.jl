using Test
import WavePropBase as WPB

@testset "Disk" begin
    disk = WPB.Disk() # abstract entity
    Γ    = WPB.boundary(disk) |> WPB.Domain
    M    = WPB.meshgen(Γ,(10,))
    M    = WPB.meshgen(Γ;meshsize=0.1)
    # plot(M,Γ)
    @test WPB.entities(Γ) == WPB.boundary(disk)
    @test WPB.geometric_dimension(disk) == 2
end

@testset "Box" begin
    box = WPB.Box() # abstract entity
    Γ    = WPB.boundary(box) |> WPB.Domain
    M    = WPB.meshgen(Γ,(10,10))
    M    = WPB.meshgen(Γ;meshsize=0.1)
    # plot(M,Γ)
    @test WPB.entities(Γ) == WPB.boundary(box)
    @test WPB.geometric_dimension(box) == 3
end

@testset "Ball" begin
    ball = WPB.Ball() # abstract entity
    Γ = WPB.boundary(ball) |> WPB.Domain
    M = WPB.meshgen(Γ,(1,1))
    M = WPB.meshgen(Γ;meshsize=0.1)
    # plot(M,Γ)
    @test WPB.entities(Γ) == WPB.boundary(ball)
    @test WPB.geometric_dimension(ball) == 3
end
