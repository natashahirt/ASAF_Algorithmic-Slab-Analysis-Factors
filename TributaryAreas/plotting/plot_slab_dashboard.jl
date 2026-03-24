"""
    plot_forces(results::SlabOptimResults, map_variable::Symbol; highlight=-1, name::String="", plot_slab_context::Bool=false, override::Bool=false)

Plots the forces on a slab based on the optimization results. The function can highlight specific elements and optionally plot the slab context.

# Arguments
- `results::SlabOptimResults`: The results from the slab optimization.
- `map_variable::Symbol`: The variable to map (e.g., :A, :displacement).
- `highlight`: Index of the element to highlight. Default is -1 (no highlight).
- `name::String`: Optional name for the plot.
- `plot_slab_context::Bool`: Whether to plot the slab context.
- `override::Bool`: Whether to override the highlight behavior.
"""
function plot_forces(results::SlabOptimResults, map_variable::Symbol; highlight=-1, name::String="", plot_slab_context::Bool=false, override::Bool=false)
    self = results.self
    model = self.model

    fig = determine_figure_size(highlight, override)
    colormap = select_colormap(map_variable)
    fontsize, smallfontsize = 11, 8

    x_positions, y_positions = extract_positions(model)
    x_range, y_range, interval = calculate_ranges(x_positions, y_positions)

    if highlight <= 0 && !override
        ax_map = setup_main_axis(fig, x_positions, y_positions, x_range, y_range, interval, fontsize)
        if plot_slab_context
            self.plot_analysis = true
            #slab = analyze_slab(self, axis=ax_map)
        end
        ax_map, limits = plot_map!(ax_map, results, map_variable, colormap, highlight=-1)
    else
        highlight = adjust_highlight(highlight, override)

        grid = GridLayout(fig[1,1])
        mapgrid = GridLayout(grid[1,1])
        subgrid = GridLayout(grid[2,1])
    
        rowsize!(grid, 1, Fixed(332.5))
        rowgap!(grid, 1, Fixed(40))
    
        ax_map = setup_main_axis(mapgrid, x_positions, y_positions, x_range, y_range, interval, fontsize)
        ax_map, limits = plot_map!(ax_map, results, map_variable, colormap, highlight=highlight, override=true)
        ax_displacement = Axis(subgrid[1,1], xlabel="x [m]", ylabel="Δ [m]", title="Global Displacement", yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, topspinevisible=false, rightspinevisible=false, titlesize=fontsize)
        ax_moment = Axis(subgrid[1,2], xlabel="x [m]", ylabel="M [Nm]", title="Moment", yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, topspinevisible=false, rightspinevisible=false, titlesize=fontsize)
        ax_shear = Axis(subgrid[1,3], xlabel="x [m]", ylabel="V [kN]", title="Shear", yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, topspinevisible=false, rightspinevisible=false, titlesize=fontsize)
        ax_load = Axis(subgrid[1,4], xlabel="x [m]", ylabel="Load [kN]", title="Load Profile (z)", yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, topspinevisible=false, rightspinevisible=false, titlesize=fontsize)
    
        linkxaxes!(ax_displacement, ax_moment, ax_shear, ax_load)
    
        lines!(ax_displacement, results.x[highlight], results.Δ_global[highlight])
        lines!(ax_moment, results.x[highlight], results.My[highlight])
        lines!(ax_shear, results.x[highlight], results.Vy[highlight])
        barplot!(ax_load, results.Px[1], results.P[1], color=results.P[1], colormap=:blues)
    end

    finalize_plot(ax_map, self.slab_type, self.vector_1d)
    return fig
end

function determine_figure_size(highlight, override)
    return if highlight <= 0 && !override
        Figure(size=(900,900))
    else
        Figure(size=(760,522.5))
    end
end

function select_colormap(map_variable)
    return map_variable == :A ? :blues : Reverse(:blues)
end

function extract_positions(model)
    x_positions = [node.position[1] for node in model.nodes]
    y_positions = [node.position[2] for node in model.nodes]
    return x_positions, y_positions
end

function calculate_ranges(x_positions, y_positions)
    x_range = maximum(x_positions) - minimum(x_positions)
    y_range = maximum(y_positions) - minimum(y_positions)
    interval = 2
    return x_range, y_range, interval
end

function setup_main_axis(fig, x_positions, y_positions, x_range, y_range, interval, fontsize)
    return Axis(fig[1,1], xlabel="x [m]", ylabel="y [m]", aspect=DataAspect(),
                xticks=(minimum(x_positions):interval:maximum(x_positions), [string(Int(i)) for i in 0:interval:x_range]),
                yticks=(minimum(y_positions):interval:maximum(y_positions), [string(Int(i)) for i in 0:interval:y_range]),
                xgridvisible=false, ygridvisible=false, titlesize=fontsize, topspinevisible=false, leftspinevisible=false,
                bottomspinevisible=false, rightspinevisible=false)
end

function adjust_highlight(highlight, override)
    return if highlight < 0 && override
        abs(highlight)
    else
        highlight
    end
end


function finalize_plot(ax_map, slab_type, vector_1d)
    title = slab_type == :isotropic ? string(slab_type) : string(slab_type) * " " * string(vector_1d)
    hidedecorations!(ax_map)
    ax_map.bottomspinevisible = false
    ax_map.leftspinevisible = false
    ax_map.alignmode = Inside()
end

"""
    plot_map!(axis::Axis, results::SlabOptimResults, variable::Symbol, colormap; highlight=-1, override::Bool=false)

Plots the map of a specified variable on the given axis.

# Arguments
- `axis::Axis`: The axis to plot on.
- `results::SlabOptimResults`: The results from the slab optimization.
- `variable::Symbol`: The variable to plot (e.g., :displacement, :moment).
- `colormap`: The colormap to use for plotting.
- `highlight`: Index of the element to highlight. Default is -1 (no highlight).
- `override::Bool`: Whether to override the highlight behavior.
"""
function plot_map!(axis::Axis, results::SlabOptimResults, variable::Symbol, colormap; highlight=-1, override::Bool=false)
    self = results.self
    model = self.model
    elements = model.elements[:beam]

    if override
        highlight *= -1
    end

    variable_min, variable_max = determine_variable_range(results, variable)

    for (i, element) in enumerate(elements)
        linewidth = determine_linewidth(highlight, i, results.areas[i])
        plotting_variable = select_plotting_variable(variable, results, i, axis)

        resolution = length(plotting_variable)
        x_series, y_series, pos = interpolate_between_nodes(element, resolution)

        x_series_doubled, y_series_doubled = double_up_series(x_series, y_series, element)

        rotation = calculate_rotation(element)
        linesegments!(axis, x_series_doubled, y_series_doubled, linewidth=linewidth, color=plotting_variable, colorrange=(variable_min, variable_max), colormap=colormap, alpha=0.3)
    end

    return axis, (variable_min, variable_max)
end

function determine_variable_range(results, variable)
    if variable in [:displacement, :Δ, :D]
        return (minimum([minimum(d) for d in results.Δ_global]), maximum([maximum(abs.(d)) for d in results.Δ_global]))
    elseif variable in [:moment, :M]
        return (minimum([minimum(m) for m in results.My]), maximum([maximum(abs.(m)) for m in results.My]))
    elseif variable in [:shear, :V]
        return (minimum([minimum(s) for s in results.Vy]), maximum([maximum(abs.(s)) for s in results.Vy]))
    elseif variable == :A
        return (minimum(results.areas), maximum(results.areas))
    else
        return (0, 0)
    end
end

function determine_linewidth(highlight, i, area)
    return if highlight <= 0
        sqrt(area)
    else
        highlight == i ? 5 : 1
    end
end

function select_plotting_variable(variable, results, i, axis)
    if variable in [:displacement, :Δ, :D]
        axis.title = "Global Displacement"
        return results.Δ_global[i]
    elseif variable in [:moment, :M]
        axis.title = "Moment"
        return results.My[i]
    elseif variable in [:shear, :V]
        axis.title = "Shear"
        return results.Vy[i]
    elseif variable == :A
        axis.title = "Sections"
        return results.areas[i]
    end
end

function double_up_series(x_series, y_series, element)
    x_series_doubled = [element.nodeStart.position[1]]
    y_series_doubled = [element.nodeStart.position[2]]

    for j in 1:lastindex(x_series)-1
        append!(x_series_doubled, [x_series[j], x_series[j]])
        append!(y_series_doubled, [y_series[j], y_series[j]])
    end

    push!(x_series_doubled, element.nodeEnd.position[1])
    push!(y_series_doubled, element.nodeEnd.position[2])

    return x_series_doubled, y_series_doubled
end

function calculate_rotation(element)
    rotation = mod(angle_pointpoint(element.nodeEnd, element.nodeStart, degrees=false, clockwise=false), pi)
    if rotation > pi/2
        rotation = pi + rotation
    end
    return rotation
end