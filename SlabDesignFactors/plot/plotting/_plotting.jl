using DataFrames
using CSV
using Colors
using CairoMakie
using Statistics
using PrettyTables

include("../../core/_core.jl")
include("assemble_data.jl")
include("clean_data.jl")
include("utils.jl")

# plot files
include("0_max_deflection.jl")
include("0_plot_slab.jl")
include("1_multiplot.jl")
include("2_megaplot.jl")
include("3_topology.jl")
include("4_surface.jl")
include("5_beam_sizes.jl")
include("6_depth.jl")
include("7_fix_params.jl")
include("8_stats_summary.jl")
include("9_stats_topology.jl")
include("10_subplots.jl")
include("11_geometries.jl")
include("12_constrained_inventory.jl")