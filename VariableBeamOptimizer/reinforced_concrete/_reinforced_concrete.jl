using Pkg
Pkg.activate(".")

using Revise

# Load necessary packages for structural analysis and optimization
using Asap, AsapToolkit, AsapOptim
using CairoMakie
using Nonconvex, Zygote
using JuMP: JuMP
using GLPK, Ipopt, Gurobi
using Statistics, Colors, DataFrames, CSV, JSON, Interpolations, StatsBase, UnPack, Dates, XLSX, Roots

# Load packages for websocket
using HTTP

# Load optimization methods
Nonconvex.@load MMA
Nonconvex.@load NLopt
Nonconvex.@load Juniper

const gurobi_env = Gurobi.Env()  # Create once and reuse

# Include necessary modules for slab analysis and optimization
include("../../SlabDesignFactors/core/_core.jl")
include("../../SlabDesignFactors/plot/plotting/_plotting.jl")
include("../../TributaryAreas/TributaryAreas.jl")
include("../../VariableBeamOptimizer/VariableBeamOptimizer.jl")

Pkg.status()

include("utils.jl")
include("utils_shear.jl")
include("utils_tension.jl")
include("evaluation.jl")
include("plotting.jl")
include("sizing.jl")
include("sizing_zonal.jl")
include("rebar_placement.jl")
include("optimization.jl")
include("../setup/_setup.jl")