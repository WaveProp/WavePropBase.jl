using Test
import WavePropBase as WPB
using StaticArrays

@testset "ReferenceLine" begin
    d = WPB.ReferenceLine()
    k = 3
    sp = WPB.PolynomialSpace(d, k)
    @test WPB.dimension(sp) == k + 1
    b = WPB.monomial_basis(sp)
    @test length(b) == WPB.dimension(sp)

    # P0 basis
    sp = WPB.PolynomialSpace(d, 0)
    nodes = [0.5]
    b = WPB.lagrange_basis(nodes, sp)
    @assert length(b) == length(nodes)
    @test b[1](0.5) ≈ 1

    # P1 basis
    sp = WPB.PolynomialSpace(d, 1)
    nodes = [0, 1]
    b = WPB.lagrange_basis(nodes, sp)
    @assert length(b) == length(nodes)
    @test b[1](nodes[1]) ≈ 1
    @test b[1](nodes[2]) < 1e-15
    @test b[2](nodes[1]) < 1e-15
    @test b[2](nodes[2]) ≈ 1
end

@testset "ReferenceTriangle" begin
    d = WPB.ReferenceTriangle()
    k = 3
    sp = WPB.PolynomialSpace(d, k)
    @test WPB.dimension(sp) == (k + 1) * (k + 2) / 2
    b = WPB.monomial_basis(sp)
    @test length(b) == WPB.dimension(sp)

    # P0 basis over triangle
    sp = WPB.PolynomialSpace(d, 0)
    nodes = [SVector(1 / 3, 1 / 3)]
    b = WPB.lagrange_basis(nodes, sp)
    @assert length(b) == length(nodes)
    @test b[1](nodes[1]) ≈ 1

    # P1 basis over triangle
    sp = WPB.PolynomialSpace(d, 1)
    nodes = [SVector(0, 0), SVector(0, 1), SVector(1, 0)]
    b = WPB.lagrange_basis(nodes, sp)
    @assert length(b) == length(nodes)
    for i in 1:length(nodes)
        for j in 1:length(nodes)
            if i == j
                @test b[i](nodes[j]) ≈ 1
            else
                @test abs(b[i](nodes[j])) < 1e-10
            end
        end
    end
end

@testset "ReferenceSquare" begin
    d = WPB.ReferenceSquare()
    k = 3
    sp = WPB.PolynomialSpace(d, k)
    @test WPB.dimension(sp) == (k + 1)^2
    b = WPB.monomial_basis(sp)
    @test length(b) == WPB.dimension(sp)
end
