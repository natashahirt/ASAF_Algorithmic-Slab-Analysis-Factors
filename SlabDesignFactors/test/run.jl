println("Loading project via _scripts.jl...")
include("../scripts/_scripts.jl")

println("Project loaded. Running tests...")
include("runtests.jl")
println("All tests complete.")
