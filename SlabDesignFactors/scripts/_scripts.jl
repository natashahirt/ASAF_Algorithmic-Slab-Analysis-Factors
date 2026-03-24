using Pkg
Pkg.activate(".")

include(joinpath(@__DIR__, "..", "makie_skip_env.jl"))
const _SLAB_SKIP_MAKIE = _slab_sync_skip_makie_env!()

if _SLAB_SKIP_MAKIE
    @info "SLABDESIGN_SKIP_MAKIE: skipping CairoMakie and SlabDesignFactors plotting (headless cluster)."
end

using Revise

# Load necessary packages for structural analysis and optimization
using Asap, AsapToolkit, AsapOptim
if !_SLAB_SKIP_MAKIE
    using CairoMakie
end
using Nonconvex, Zygote
using JuMP: JuMP
using GLPK, Ipopt, Gurobi, NLopt
using Statistics, Colors, DataFrames, CSV, JSON, Interpolations, StatsBase, UnPack, Dates, XLSX

# Load packages for websocket
using HTTP

# Load optimization methods
Nonconvex.@load MMA
Nonconvex.@load NLopt
Nonconvex.@load Juniper

# Include necessary modules for slab analysis and optimization
include("../../SlabDesignFactors/core/_core.jl")
if !_SLAB_SKIP_MAKIE
    include("../../SlabDesignFactors/plot/plotting/_plotting.jl")
end
include("../../TributaryAreas/TributaryAreas.jl")
include("../../VariableBeamOptimizer/VariableBeamOptimizer.jl")
# Projects/ (IndeterminateTopopt, etc.) pulls in CairoMakie-heavy code; omitted so tests and headless
# runs work. From a REPL with the same working directory as this script: `include("../../Projects/_projects.jl")`.

Pkg.status()
