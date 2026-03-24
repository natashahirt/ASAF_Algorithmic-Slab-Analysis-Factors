using Pkg
Pkg.activate(".")

"""
Headless / no-Makie mode for cluster login nodes (Linux without DISPLAY/WAYLAND).
Set `SLABDESIGN_SKIP_MAKIE=1` to force, or `=0` to force CairoMakie (e.g. X11 forwarding).
"""
function _slab_skip_makie()
    if haskey(ENV, "SLABDESIGN_SKIP_MAKIE")
        v = strip(get(ENV, "SLABDESIGN_SKIP_MAKIE", ""))
        isempty(v) && return _slab_infer_headless_skip_makie()
        lv = lowercase(v)
        if lv in ("0", "false", "no", "n", "off")
            return false
        end
        return lv in ("1", "true", "yes", "y", "on")
    end
    return _slab_infer_headless_skip_makie()
end

function _slab_infer_headless_skip_makie()
    Sys.iswindows() && return false
    Sys.isapple() && return false
    return isempty(get(ENV, "DISPLAY", "")) && isempty(get(ENV, "WAYLAND_DISPLAY", ""))
end

const _SLAB_SKIP_MAKIE = _slab_skip_makie()
ENV["SLABDESIGN_SKIP_MAKIE"] = _SLAB_SKIP_MAKIE ? "1" : "0"

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
