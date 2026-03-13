using Pkg
Pkg.activate(".")

using Revise

# Load necessary packages for structural analysis and optimization
using Asap, AsapToolkit, AsapOptim
using CairoMakie
using Nonconvex, Zygote
using JuMP: JuMP
using GLPK, Ipopt, Gurobi
using Statistics, Colors, DataFrames, CSV, JSON, Interpolations, StatsBase, UnPack, Dates, XLSX

# Load packages for websocket
using HTTP

# Load optimization methods
Nonconvex.@load MMA
Nonconvex.@load NLopt
Nonconvex.@load Juniper

# Include necessary modules for slab analysis and optimization
include("../../SlabDesignFactors/core/_core.jl")
include("../../SlabDesignFactors/plot/plotting/_plotting.jl")
include("../../TributaryAreas/TributaryAreas.jl")
include("../../VariableBeamOptimizer/VariableBeamOptimizer.jl")
include("../../Projects/_projects.jl")

Pkg.status()