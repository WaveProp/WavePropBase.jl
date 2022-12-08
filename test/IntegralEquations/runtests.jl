using SafeTestsets

@safetestset "NystromMesh" begin
    include("nystrommesh_test.jl")
end

@safetestset "NystromMesh" begin
    include("dim_test.jl")
end

@safetestset "Adaptive integration with HCubature" begin
    include("hcubature_test.jl")
end
