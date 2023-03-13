using SafeTestsets

@safetestset "NystromMesh" begin
    include("nystrommesh_test.jl")
end

@safetestset "NystromMesh" begin
    include("dim_test.jl")
end

@safetestset "Adaptive integration with HCubature" begin
    # TODO: include this tests once a new version of HCubature with my PR gets released
    # include("hcubature_test.jl")
end
