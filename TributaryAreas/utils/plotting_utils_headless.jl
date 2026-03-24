# Stubs when CairoMakie is not loaded (e.g. SLABDESIGN_SKIP_MAKIE=1 on headless nodes).
# Analysis code with plot_analysis=false never calls these; they exist only for load order / clear errors.

function _tributary_plotting_disabled!()
    error(
        "TributaryAreas plotting requires CairoMakie. " *
        "Unset SLABDESIGN_SKIP_MAKIE or set it to 0, or run with plot=false / plot_analysis=false.",
    )
end

function setup_plot(plot::Bool=true; size=(1200, 800))
    _tributary_plotting_disabled!()
end

function plot_trimmed_line(ax, trimmed_line, trimmed_distance; color=nothing)
    _tributary_plotting_disabled!()
end

function plot_intersection(self::SlabAnalysisParams, x, y, point, distance; color=nothing)
    _tributary_plotting_disabled!()
end

function plot_node!(ax, args...; kwargs...)
    _tributary_plotting_disabled!()
end

function plot_elements!(ax, elements; kwargs...)
    _tributary_plotting_disabled!()
end

function plot_line!(ax, args...; kwargs...)
    _tributary_plotting_disabled!()
end

function plot_ray!(ax, node::Vector{Float64}, direction::Vector{Float64}, magnitude::Float64; color=:black)
    _tributary_plotting_disabled!()
end

function plot_model(model::Asap.Model; numbers::Bool=false, plot_context::PlotContext, ax=nothing, type_information::Dict=Dict())
    _tributary_plotting_disabled!()
end

function plot_model(model::Asap.Model)
    _tributary_plotting_disabled!()
end
