"""
    find_opposite_edge_orthogonal(point::Vector{Float64}, edge_list::Vector; vector_1d=[1.0, 0.0], plot=true, color::Union{Nothing, Symbol}=nothing)

Finds the opposite edge orthogonal to a given point and vector direction.

# Arguments
- `point::Vector{Float64}`: The starting point for the ray.
- `edge_list::Vector`: List of edges to check for intersection.
- `vector_1d`: Direction vector for the ray (default is [1.0, 0.0]).
- `plot`: Boolean to determine if the intersection should be plotted (default is true).
- `color::Union{Nothing, Symbol}`: Color for plotting (default is nothing).

# Returns
- `target_point`: The intersection point on the target edge.
- `distance`: The distance from the point to the intersection.
"""
function find_opposite_edge_orthogonal(self::SlabAnalysisParams, point::Vector{Float64}, edge_list::Vector; vector_1d=self.vector_1d, color=nothing, plot=true)
    for j in 1:lastindex(edge_list)
        target_edge = edge_list[j]
        target_point, distance = ray_line_intersect(point, vector_1d, target_edge, param=1.0)

        if distance == 0
            continue
        end

        opacity = self.slab_type == :orth_biaxial ? 0.5 : 1.0

        if self.plot_context.plot && plot
            max_distance = 3
            normalized_param = distance / max_distance
            plot_line!(self.plot_context.ax, [point[1], target_point[1]], [point[2], target_point[2]], color=(isnothing(color) ? normalized_param : color), alpha=opacity)
        end

        return target_point, distance
    end

    return [0.0, 0.0], 0.0

end

"""
    find_opposite_edge(point::Vector{Float64}, i::Int64, edge_list::Vector{Vector}, element_list::Vector{T}; vector_1d=[1.0, 0.0], fix_param=false, plot=false, orthogonal=false, model_elements::Vector{Element{T}}) where T

Finds the opposite edge for a given point and element index.

# Arguments
- `point::Vector{Float64}`: The starting point for the ray.
- `i::Int64`: Index of the current element.
- `edge_list::Vector{Vector}`: List of edges to check for intersection.
- `element_list::Vector{T}`: List of elements corresponding to edges.
- `vector_1d`: Direction vector for the ray (default is [1.0, 0.0]).
- `fix_param`: Boolean to fix the parameter for stiffness ratio (default is false).
- `plot`: Boolean to determine if the intersection should be plotted (default is false).
- `orthogonal`: Boolean to determine if orthogonal intersection is required (default is false).
- `model_elements::Vector{Element{T}}`: List of model elements.

# Returns
- `distance`: The distance to the intersection.
- `target_point`: The intersection point on the target edge (if orthogonal).
- `target_element`: The target element (if orthogonal).
"""
function find_opposite_edge(self::SlabAnalysisParams, point::Vector{Float64}, i::Int64, edge_list::Vector{Vector}, element_list; orthogonal=false, model_elements)
    single_element = i <= 0

    for j in 1:lastindex(edge_list)
        if i == j
            continue # Skip if the same element
        end

        target_edge = edge_list[j]
        target_element = element_list[j]

        stiffness_ratio = if single_element
            1.0
        elseif self.fix_param
            0.5
        else
            get_stiffness_ratio(element_list[i], target_element)
        end

        vector_1d = self.perp == true ? self.perp_vector_1d : self.vector_1d

        target_point, distance = ray_line_intersect(point, vector_1d, target_edge, param=stiffness_ratio)

        if distance == 0
            continue
        end

        if self.plot_context.plot
            max_distance = 3
            normalized_param = distance / max_distance
            opacity = self.slab_type == :orth_biaxial ? 0.5 : 1.0
            plot_line!(self.plot_context.ax, [point[1], target_point[1]], [point[2], target_point[2]], color=normalized_param, alpha=opacity)
        end

        if orthogonal
            return target_point, distance, target_element
        end

        return distance
    end
end

# STRAIGHT SKELETON

"""
    get_next_wavefront(bisectors::Vector{Vector}, startpoints::Vector{Vector}; distance::Float64=0.1)

Calculates the next wavefront points based on given bisectors and startpoints.

# Arguments
- `bisectors::Vector{Vector}`: A vector of bisector vectors.
- `startpoints::Vector{Vector}`: A vector of starting points corresponding to each bisector.
- `distance::Float64`: The distance to move along each bisector (default is 0.1).

# Returns
- `bisectors`: The original bisectors.
- `new_points`: The new points calculated by moving along each bisector from the startpoints.
"""
function get_next_wavefront(bisectors::Vector{Vector}, startpoints::Vector{Vector}; distance::Float64=0.1)
    new_points = Vector{Vector}()

    for i in 1:lastindex(bisectors)
        if bisectors[i] == [0.0, 0.0]
            # If the bisector is zero, retain the original startpoint
            push!(new_points, startpoints[i])
        else
            # Scale the bisector and calculate the new point
            scaled_bisector = unit_vector(bisectors[i]) * distance
            point = startpoints[i] + scaled_bisector
            push!(new_points, point)
        end
    end

    return bisectors, new_points
end
