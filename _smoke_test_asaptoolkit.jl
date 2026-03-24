using Pkg
const ROOT = "C:/Users/nhirt/MIT Dropbox/Natasha Hirt/_MIT/DS Research/2024_Slab Design Factors"
Pkg.activate(ROOT)
Pkg.instantiate()
println("=== instantiate OK ===")

using Asap, AsapToolkit, AsapOptim
println("=== AsapToolkit loaded ===")

include(joinpath(ROOT, "TributaryAreas", "TributaryAreas.jl"))
println("=== TributaryAreas loaded ===")

include(joinpath(ROOT, "VariableBeamOptimizer", "VariableBeamOptimizer.jl"))
println("=== VariableBeamOptimizer loaded ===")

include(joinpath(ROOT, "SlabDesignFactors", "core", "_core.jl"))
println("=== SlabDesignFactors core loaded ===")

println("ALL SMOKE TESTS PASSED")
