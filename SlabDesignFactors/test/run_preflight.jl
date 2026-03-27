include(joinpath(@__DIR__, "_license_env.jl"))

using Test

println("Loading project via _scripts.jl for preflight tests...")
include("../scripts/_scripts.jl")

println("Project loaded. Running preflight parity tests...")

@testset "SlabDesignFactors Preflight" begin
    include("test_preflight_mip_nlp.jl")
end

println("Preflight tests complete.")
