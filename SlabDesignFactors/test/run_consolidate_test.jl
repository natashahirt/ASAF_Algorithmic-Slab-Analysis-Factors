include(joinpath(@__DIR__, "_license_env.jl"))

println("Loading project...")
include("../scripts/_scripts.jl")

println("Running consolidate_loads equilibrium test...")
include("test_consolidate_loads.jl")
println("Done.")
