module SlabDesignFactors

using Pkg
Pkg.activate(".")

include(joinpath(@__DIR__, "makie_skip_env.jl"))
const _SLAB_SKIP_MAKIE = _slab_sync_skip_makie_env!()
if _SLAB_SKIP_MAKIE
    @info "SLABDESIGN_SKIP_MAKIE: skipping CairoMakie (headless, forced env, or cluster default)."
end

using Revise

# Load necessary packages for structural analysis and optimization
using Asap, AsapToolkit, AsapOptim
if !_SLAB_SKIP_MAKIE
    using CairoMakie
end
using Nonconvex, Zygote
using JuMP: JuMP
using GLPK, Ipopt, NLopt, Gurobi
using Statistics, Colors, DataFrames, CSV, JSON, Interpolations, StatsBase, UnPack

# Load packages for websocket
using HTTP
# Load optimization methods
Nonconvex.@load MMA
Nonconvex.@load NLopt

# Include necessary modules for slab analysis and optimization
include("../SlabDesignFactors/core/_core.jl")
include("../TributaryAreas/TributaryAreas.jl")
include("../VariableBeamOptimizer/VariableBeamOptimizer.jl")

Pkg.status()

println("Initialization of Slab_Design_Factors complete")

end # module SlabDesignFactors
