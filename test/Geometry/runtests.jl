using Test
using SafeTestsets

@safetestset "Reference shapes tests" begin
    include("referenceshapes_test.jl")
end

@safetestset "Domain tests" begin
    include("domain_test.jl")
end
