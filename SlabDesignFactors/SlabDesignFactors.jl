module SlabDesignFactors

using Pkg
Pkg.activate(".")

using Revise

# Load necessary packages for structural analysis and optimization
using Asap, AsapToolkit, AsapOptim
using CairoMakie
using Nonconvex, Zygote
using JuMP: JuMP
using GLPK, Ipopt
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
