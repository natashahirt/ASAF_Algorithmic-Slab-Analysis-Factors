"""
    setup_plot(axis)

Sets up a plot for geometry generation. If no axis is provided, it creates a new figure and axis with a specified size and aspect ratio.

# Arguments
- `ax::Union{Axis, Nothing}`: An existing axis to use for plotting. If `nothing`, a new axis is created.

# Returns
- `ax`: The axis to be used for plotting.

# Usage
This function is typically used in the context of visualizing geometric configurations or results in a graphical format.
"""
function setup_plot(plot::Bool = true; size=(1200, 800))
    fig = Figure(size = size);
    ax = Axis(fig[1, 1], aspect = DataAspect());
    plot_context = PlotContext(plot, fig, ax)

    Colorbar(plot_context.fig[1, 2], halign=:right, tellheight = false, height=Relative(0.5), colormap=:Blues_9, limits=(0,3), label = "Distance [m]", labelsize=20, highclip=:black)
    hidespines!(plot_context.ax)
    plot_context.ax.xticklabelsize = 20
    plot_context.ax.yticklabelsize = 20
    return plot_context
end


"""
    plot_trimmed_line(trimmed_line, trimmed_distance, color)

Plots a line segment with a specified color, representing a trimmed section of a line.

# Arguments
- `ax::Axis`: The axis on which to plot.
- `trimmed_line::Vector{Vector{Float64}}`: The endpoints of the line segment to plot.
- `trimmed_distance::Float64`: The length of the trimmed line segment.
- `color::Union{Nothing, Symbol}`: The color to use for plotting. If `nothing`, a normalized parameter based on distance is used.

# Usage
This function is used to visualize trimmed line segments, typically in the context of geometric configurations.
"""
function plot_trimmed_line(ax::Axis, trimmed_line, trimmed_distance; color=nothing)
    max_distance = 3
    normalized_param = trimmed_distance / max_distance
    if !isnothing(color)
        plot_line!(ax, [trimmed_line[1][1], trimmed_line[2][1]], [trimmed_line[1][2], trimmed_line[2][2]], color=color)
    else
        plot_line!(ax, [trimmed_line[1][1], trimmed_line[2][1]], [trimmed_line[1][2], trimmed_line[2][2]], color=normalized_param)
    end
end

"""
    plot_intersection(ax, x, y, point, distance, color)

Plots a line from a given point to an intersection point with a specified color.

# Arguments
- `ax`: The axis on which to plot.
- `x::Float64`: The x-coordinate of the starting point.
- `y::Float64`: The y-coordinate of the starting point.
- `point::Vector{Float64}`: The intersection point as a vector of coordinates.
- `distance::Float64`: The distance from the starting point to the intersection point.
- `color::Union{Nothing, Symbol}`: The color to use for plotting. If `nothing`, a normalized parameter based on distance is used.

# Usage
This function is used to visualize the intersection of lines or paths in a plot, typically in the context of geometric configurations.
"""
function plot_intersection(self::SlabAnalysisParams, x, y, point, distance; color=nothing)
    max_distance = 3
    normalized_param = distance / max_distance
    if !isnothing(color)
        plot_line!(self.plot_context.ax, [x, point[1]], [y, point[2]], color=color)
    else
        plot_line!(self.plot_context.ax, [x, point[1]], [y, point[2]], color=normalized_param)
    end
end



"""
    plot_node!(ax::Axis, points::Vector{Vector}; color=:black, marker=:circle, numbers=false)

Plot nodes on the given axis `ax` using the specified `color` and `marker`.
Optionally, enumerate the nodes if `numbers` is true.
"""
function plot_node!(ax::Axis, points::Vector{Vector}; color=:black, marker=:circle, numbers=false)
    for (i, point) in enumerate(points)
        plot_node!(ax, point, color=color, marker=marker)
        if numbers
            label = string(i)
            text!(ax, label, position=(point[1], point[2] + 0.1), align=(:center, :bottom), color=:black)
        end
    end
end

"""
    plot_node!(ax::Axis, point::Vector{Float64}; color=:black, marker=:circle)

Plot a single node on the given axis `ax` using the specified `color` and `marker`.
"""
function plot_node!(ax::Axis, point::Vector{Float64}; color=:black, marker=:circle)
    plot_node!(ax, point[1], point[2], color=color, marker=marker)
end

"""
    plot_node!(ax::Axis, x::Float64, y::Float64; color=:black, marker=:circle)

Plot a single node at coordinates `(x, y)` on the given axis `ax` using the specified `color` and `marker`.
"""
function plot_node!(ax::Axis, x::Float64, y::Float64; color=:black, marker=:circle)
    scatter!(ax, x, y, color=color, marker=marker, strokewidth=1, strokecolor=:white)
end

"""
    plot_node!(ax::Axis, node::Node; color=:black)

Plot a node on the given axis `ax` with properties determined by the node's degrees of freedom (`dof`) and `id`.
"""
function plot_node!(ax::Axis, node::Node; color=:black)
    x, y = node.position[1], node.position[2]
    marker = :circle

    if node.dof == [0, 0, 0, 0, 0, 0] || node.id == :fixed
        return
    elseif node.dof == [0, 0, 0, 1, 1, 1]
        marker = :utriangle
    elseif node.dof == [1, 1, 1, 0, 0, 0] || node.dof == [1, 1, 1, 1, 1, 1]
        marker = :cross
        color = :deeppink
    end

    if node.id == :column
        marker = :circle
        color = :deeppink
    elseif node.id == :wall
        marker = :circle
        color = :black
    end

    plot_node!(ax, x, y, color=color, marker=marker)
end

"""
    plot_node!(ax::Axis, nodes::Vector{Node}; color=:black, numbers=false)

Plot multiple nodes on the given axis `ax` using the specified `color`.
Optionally, enumerate the nodes if `enumerate` is true.
"""
function plot_node!(ax::Axis, nodes::Vector{Node}; color=:black, numbers=false)
    for (i, node) in enumerate(nodes)
        plot_node!(ax, node, color=color)
        if numbers 
            label = node.id == :fixed ? "" : string(node.nodeID)
            text!(ax, label, position=(node.position[1], node.position[2] + 0.1), align=(:center, :bottom), color=:black)
        end
    end
end

"""
    plot_elements!(ax::Axis, elements::Vector{<:Element}; color=:black)

Plot elements on the given axis `ax` using the specified `color`.
"""
function plot_elements!(ax::Axis, elements::Union{Vector{Element{T}}, Vector{Element}}; color=:black, linewidth=2, type_information::Dict=Dict()) where T <: Asap.Release
    for (i, element) in enumerate(elements)
        if element.nodeStart.id == :wall && element.nodeEnd.id == :wall
            linewidth = 5
            line_color = :grey
        elseif typeof(element) == Element{Asap.FreeFree}
            linewidth = 2
            line_color = :deeppink
        else
            linewidth = 2
            line_color = :black
        end

        if haskey(type_information, "i_perimeter") && i in type_information["i_perimeter"]
            linestyle = :dash   
        else
            linestyle = :solid
        end
        plot_line!(ax, [element.nodeStart, element.nodeEnd], color=line_color, linewidth=linewidth, linestyle=linestyle)
    end
end

"""
    plot_line!(ax::Axis, x_coords::Vector{Float64}, y_coords::Vector{Float64}; color=1)

Plot a line on the given axis `ax` using the specified `color`.
"""
function plot_line!(ax::Axis, x_coords::Vector{Float64}, y_coords::Vector{Float64}; color=1, alpha=1.0, linewidth=5, linestyle=:solid)
    if isa(color, Symbol) || isa(color, Color)
        lines!(ax, x_coords, y_coords, color=color, linewidth=linewidth, alpha=alpha, linestyle=linestyle)
    else
        line_color = if color < 0
            :red
        elseif color > 1
            :black
        else
            color
        end

        lines!(ax, x_coords, y_coords, color=line_color, colormap=:Blues_9, colorrange=(0, 1), linewidth=linewidth, alpha=alpha, linestyle=linestyle) # linewidth=3 for larger slabs
    end
end

"""
    plot_line!(ax::Axis, nodes::Vector{Node}; closed=false, color=:black)

Plot a line connecting nodes on the given axis `ax`.
If `closed` is true, the line will form a closed loop.
"""
function plot_line!(ax::Axis, nodes::Vector{Node}; closed=false, color=:black, linewidth=2, linestyle=:solid)
    x_coords = [node.position[1] for node in nodes]
    y_coords = [node.position[2] for node in nodes]
    if closed
        push!(x_coords, x_coords[1])
        push!(y_coords, y_coords[1])
    end
    lines!(ax, x_coords, y_coords, color=color, linewidth=linewidth, linestyle=linestyle)
end

"""
    plot_line!(ax::Axis, vector::Vector{Vector{Float64}}; color=:black)

Plot a line between two points defined by `vector` on the given axis `ax`.
"""
function plot_line!(ax::Axis, vector::Vector{Vector{Float64}}; color=:black, linewidth=2, linestyle=:solid)
    plot_line!(ax, [vector[1][1], vector[2][1]], [vector[1][2], vector[2][2]], color=color, linewidth=linewidth, linestyle=linestyle)
end

"""
    plot_ray!(ax::Axis, node::Vector{Float64}, direction::Vector{Float64}, magnitude::Float64; color=:black)

Plot a ray starting from `node` in the given `direction` with the specified `magnitude` on the axis `ax`.
"""
function plot_ray!(ax::Axis, node::Vector{Float64}, direction::Vector{Float64}, magnitude::Float64; color=:black)
    unitvector = unit_vector(direction)
    scaled_vector = unitvector * magnitude
    node_2 = node + scaled_vector
    plot_line!(ax, [node[1], node_2[1]], [node[2], node_2[2]], color=color)
end

"""
    plot_model(model::Asap.Model)

Plot the structural model.

# Arguments
- `model::Asap.Model`: The model to be plotted.
"""
function plot_model(model::Asap.Model; numbers::Bool=false, plot_context::PlotContext, ax::Union{Nothing, Axis}=nothing, type_information::Dict=Dict())
    if isnothing(ax)
        if isnothing(plot_context.ax)
            plot_context = setup_plot(plot_context.plot; size=(1200, 800))
        end
        plot_elements!(plot_context.ax, model.elements, type_information=type_information)
        plot_node!(plot_context.ax, model.nodes, numbers=numbers)
    else
        plot_elements!(ax, model.elements)
        plot_node!(ax, model.nodes, numbers=numbers)
    end
end