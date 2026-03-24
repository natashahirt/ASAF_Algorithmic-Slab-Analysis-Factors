"""
    get_interior_polygon(edges, i_nonzero_edges::Vector{Int64}, fix_param::Bool)

Calculate the interior polygon of a given set of edges. This function computes
the bisectors for each edge and determines the new edges of the polygon.

# Arguments
- `edges`: A vector of edges.
- `i_nonzero_edges`: Indices of non-zero edges.
- `fix_param`: A boolean flag to determine the method of bisector calculation.

# Returns
- `bisectors`: A vector of bisector vectors.
- `intersections`: Intersection points of bisectors.
- `new_edges`: Updated edges of the polygon.
"""
function get_interior_polygon(edges, i_nonzero_edges::Vector{Int64}, fix_param::Bool)
    bisectors = Vector[]  # Initialize an empty vector for bisectors
    startpoints = Vector[]  # Initialize an empty vector for start points
    new_edges = [edge for edge in edges]  # Copy edges to new_edges

    # Iterate over each non-zero edge index
    for i in 1:lastindex(i_nonzero_edges)
        edge, next_edge, previous_edge = get_edges(edges, i_nonzero_edges, i)
        node_a, node_o, node_b = edge[1], edge[2], next_edge[2]  # Define nodes
        push!(startpoints, node_o)  # Add node_o to startpoints

        param = calculate_param(node_a, node_o, node_b, fix_param)
        bisector_vector = calculate_bisector_vector(node_o, node_a, node_b, param)
        push!(bisectors, bisector_vector)  # Add bisector vector to bisectors
    end

    bisectors, intersections = determine_intersections(bisectors, startpoints, fix_param)    
    update_new_edges!(new_edges, i_nonzero_edges, intersections)

    return bisectors, intersections, new_edges  # Return results
end

# Subfunction to get current, next, and previous edges
function get_edges(edges, i_nonzero_edges, i)
    edge = edges[i_nonzero_edges[i]]
    next_edge = edges[i_nonzero_edges[mod1(i+1, lastindex(i_nonzero_edges))]]
    previous_edge = edges[i_nonzero_edges[mod1(i-1, lastindex(i_nonzero_edges))]]
    return edge, next_edge, previous_edge
end

# Subfunction to calculate the parameter for bisector calculation
function calculate_param(node_a, node_o, node_b, fix_param)
    if fix_param
        return 0.5
    else
        length_1 = distance_pointpoint(node_a, node_o)
        length_2 = distance_pointpoint(node_o, node_b)
        return length_1 / (length_1 + length_2)
    end
end

# Subfunction to calculate the bisector vector
function calculate_bisector_vector(node_o, node_a, node_b, param)
    bisector_point = line_line_bisector(node_o, node_a, node_b, param=param)
    if isnan(bisector_point[1]) || isnan(bisector_point[2])
        return [0.0, 0.0]  # Handle NaN values
    else
        return bisector_point - node_o  # Calculate bisector vector
    end
end

# Subfunction to determine intersections
function determine_intersections(bisectors, startpoints, fix_param)
    if fix_param
        return get_bisector_intersections(bisectors, startpoints)
    else
        return get_next_wavefront(bisectors, startpoints, distance=0.1)
    end
end

# Subfunction to update new edges with interior edges
function update_new_edges!(new_edges, i_nonzero_edges, intersections)
    interior_edges = circshift(get_cycle_edges(intersections), 1)
    for i in 1:lastindex(i_nonzero_edges)
        new_edges[i_nonzero_edges[i]] = interior_edges[i]
    end
end

"""
    get_element_params(element::Element, node_list; spacing=0.1, plot=false, orthogonal=false, color::Union{Nothing, Symbol})

Calculate parameters for an element based on its nodes and a given node list.
Interpolates points between nodes and calculates distances to the polygon.

# Arguments
- `element`: The element for which parameters are calculated.
- `node_list`: List of nodes defining the polygon.
- `spacing`: Spacing between interpolated points.
- `plot`: Boolean flag to enable plotting.
- `orthogonal`: Boolean flag to return orthogonal distances.
- `color`: Color for plotting.

# Returns
- If `orthogonal` is true: `median_points`, `element_points`, `distances`.
- Otherwise: `params`, `distances`.
"""
function get_element_params(self::SlabAnalysisParams, element::Element, node_list)
    element_length = distance_pointpoint(element.nodeStart.position[1], element.nodeEnd.position[1], element.nodeStart.position[2], element.nodeEnd.position[2])
    n = Int(floor(element_length / self.spacing))  # Determine number of interpolation points

    interp_x, interp_y, interp_params = interpolate_between_nodes(element, n)
    m_perp, c_perp = perpendicularline_equation(element.nodeStart, element.nodeEnd, element.nodeEnd)

    if abs(interp_x[1] - node_list[1][1]) > abs(interp_x[1] - node_list[end][1])
        node_list = reverse(node_list)
    end

    distances, params, median_points, element_points = calculate_distances_and_params(self, interp_x, interp_y, interp_params, node_list, m_perp)
    
    if !self.record_tributaries
        return params, distances, [0.0, 0.0]
    else
        if !isempty(element_points) && !isempty(median_points)
            vector = median_points[1] - element_points[1]
        else
            vector = [0.0, 0.0]
        end
        return params, distances, vector
    end
end

# Subfunction to calculate distances and parameters
function calculate_distances_and_params(self::SlabAnalysisParams, interp_x, interp_y, interp_params, node_list, m_perp)
    distances = Float64[]  # Initialize distances array
    params = Float64[]  # Initialize params array
    median_points = Vector[]  # Initialize median points array
    element_points = Vector[]  # Initialize element points array

    for i in 1:lastindex(interp_params)
        x, y = interp_x[i], interp_y[i]  # Current interpolated point
        intersection_points, intersection_distances = find_intersections(self, x, y, node_list, m_perp)

        if isempty(intersection_distances)
            continue
        end

        distance, point = sort_and_select_intersection(intersection_distances, intersection_points)

        if self.plot_context.plot
            plot_intersection(self, x, y, point, distance)
        end

        push!(element_points, [x, y])  # Add element point
        push!(median_points, point)  # Add median point
        push!(distances, distance)  # Add distance
        push!(params, interp_params[i])  # Add parameter
    end

    return distances, params, median_points, element_points
end

# Subfunction to find intersections
function find_intersections(self::SlabAnalysisParams, x, y, node_list, m_perp)
    intersection_points = Vector[]  # Initialize intersection points array
    intersection_distances = Float64[]  # Initialize intersection distances array

    for j in 1:lastindex(node_list) - 1

        point, distance = calculate_intersection(self, x, y, node_list[j], node_list[j+1], m_perp)
        if isnothing(point) || isnan(distance) || isinf(distance) || distance == 0
            continue
        end

        push!(intersection_points, point)  # Add intersection point
        push!(intersection_distances, distance)  # Add intersection distance
    end

    return intersection_points, intersection_distances
end

# Subfunction to calculate intersection
function calculate_intersection(self, x, y, node_start, node_end, m_perp)
    if isinf(m_perp) # flat
        m, c = line_equation(node_start, node_end)
        point_y = m * x + c
        if (point_y - minimum([node_start[2], node_end[2]])) < 0 || (maximum([node_start[2], node_end[2]]) - point_y) < 0
            return nothing, NaN
        end
        point = [x, point_y]
        distance = distance_pointpoint([x, y], point)
    elseif m_perp == 0 # vertical
        m, c = line_equation(node_start, node_end)
        if isinf(m)
            point_x = node_start[1]
        elseif isnan(m)
            point_x = node_start[1]
        else
            point_x = (y - c) / m
        end
        if !(abs(point_x - minimum([node_start[1], node_end[1]])) < 1e-6 || abs(maximum([node_start[1], node_end[1]]) - point_x) < 1e-6)
            if (point_x - minimum([node_start[1], node_end[1]])) < 0
                return nothing, NaN
            elseif (maximum([node_start[1], node_end[1]]) - point_x) < 0
                return nothing, NaN
            end
        end
        point = [point_x, y]
        distance = distance_pointpoint([x, y], point)
    else
        point, distance = ray_line_intersect([x, y], [x, m_perp * x], node_start, node_end, param=1)
    end

    return point, distance
end

# Subfunction to sort and select intersection
function sort_and_select_intersection(intersection_distances, intersection_points)
    p = sortperm(intersection_distances)
    intersection_distances, intersection_points = intersection_distances[p], intersection_points[p]
    return intersection_distances[1], intersection_points[1]
end

"""
    get_element_params_nonconvex(element::Element, node_list; spacing=0.1, plot=false, orthogonal=false, color::Union{Nothing, Symbol}=nothing)

Calculate parameters for a non-convex element based on its nodes and a given node list.
Interpolates points between nodes and calculates distances to the polygon.

# Arguments
- `ax`: The axis on which to plot.
- `element`: The element for which parameters are calculated.
- `node_list`: List of nodes defining the polygon.
- `spacing`: Spacing between interpolated points.
- `plot`: Boolean flag to enable plotting.
- `orthogonal`: Boolean flag to return orthogonal distances.
- `color`: Color for plotting.

# Returns
- `params`: Interpolated parameters.
- `distances`: Distances from interpolated points to the polygon.
"""
function get_element_params_nonconvex(self::SlabAnalysisParams, element::Element, node_list)
    # Extract element node positions
    x1_el, x2_el, y1_el, y2_el = element.nodeStart.position[1], element.nodeEnd.position[1], element.nodeStart.position[2], element.nodeEnd.position[2]
    element_vector = [x2_el - x1_el, y2_el - y1_el]  # Calculate element vector
    element_unit_vector = unit_vector(element_vector)  # Calculate unit vector

    # Project nodes onto element unit vector
    t = [dotproduct(node, element_unit_vector) for node in node_list]
    max_point, min_point = find_extreme_projection_points(node_list, t)

    # Calculate projections of element start and end points
    startpoint_t, endpoint_t = calculate_projection_points(x1_el, y1_el, x2_el, y2_el, element_unit_vector)

    # Calculate projections of max and min vectors
    min_point_proj, max_point_proj = calculate_projection_vectors(x1_el, y1_el, max_point, min_point, element_unit_vector)

    # Calculate distances
    distance_proj = calculate_distance_projection(min_point_proj, max_point_proj)

    n = Int(floor(distance_proj / self.spacing))  # Determine number of interpolation points
    interp_x, interp_y, interp_params = interpolate_between_nodes(min_point_proj, max_point_proj, n)
    m_perp, c_perp = perpendicularline_equation(element.nodeStart, element.nodeEnd, element.nodeEnd)

    # Reverse node_list if necessary
    if abs(interp_x[1] - node_list[1][1]) > abs(interp_x[1] - node_list[end][1])
        node_list = reverse(node_list)
    end

    distances, params = calculate_distances_and_params_nonconvex(self, interp_x, interp_y, interp_params, node_list, m_perp, startpoint_t, endpoint_t, element_unit_vector)

    return params, distances  # Return results
end

# Subfunction to find extreme projection points
function find_extreme_projection_points(node_list, t)
    max_point = node_list[argmax(t)]  # Find max projection point
    min_point = node_list[argmin(t)]  # Find min projection point
    return max_point, min_point
end

# Subfunction to calculate projection points
function calculate_projection_points(x1_el, y1_el, x2_el, y2_el, element_unit_vector)
    startpoint_t = dotproduct([x1_el, y1_el], element_unit_vector)
    endpoint_t = dotproduct([x2_el, y2_el], element_unit_vector)
    return startpoint_t, endpoint_t
end

# Subfunction to calculate projection vectors
function calculate_projection_vectors(x1_el, y1_el, max_point, min_point, element_unit_vector)
    max_vector = [max_point[1] - x1_el, max_point[2] - y1_el]
    min_vector = [min_point[1] - x1_el, min_point[2] - y1_el]
    min_point_proj = dotproduct(min_vector, element_unit_vector) * element_unit_vector + [x1_el, y1_el]
    max_point_proj = dotproduct(max_vector, element_unit_vector) * element_unit_vector + [x1_el, y1_el]
    return min_point_proj, max_point_proj
end

# Subfunction to calculate distance projection
function calculate_distance_projection(min_point_proj, max_point_proj)
    return distance_pointpoint(min_point_proj, max_point_proj)
end

# Subfunction to calculate distances and parameters for non-convex elements
function calculate_distances_and_params_nonconvex(self::SlabAnalysisParams, interp_x, interp_y, interp_params, node_list, m_perp, startpoint_t, endpoint_t, element_unit_vector)
    distances = Float64[]  # Initialize distances array
    params = Float64[]  # Initialize params array

    for i in 1:lastindex(interp_params)
        x, y = interp_x[i], interp_y[i]  # Current interpolated point
        projected_t = dotproduct([x, y], element_unit_vector)
        on_element = minimum([startpoint_t, endpoint_t]) <= projected_t <= maximum([startpoint_t, endpoint_t])

        intersection_points, intersection_distances = find_intersections(self, x, y, node_list, m_perp)

        if isempty(intersection_distances)
            continue
        end

        # Sort intersection distances and points
        p = sortperm(intersection_distances)
        intersection_distances, intersection_points = intersection_distances[p], intersection_points[p]

        # Check if point is on element
        if on_element
            pushfirst!(intersection_distances, 0)
            pushfirst!(intersection_points, [x, y])
        end

        total_distance = calculate_total_distance(intersection_points, intersection_distances, self.plot_context)
        push!(distances, total_distance)  # Add total distance
        push!(params, interp_params[i])  # Add parameter
    end

    return distances, params
end

# Subfunction to calculate total distance and plot if necessary
function calculate_total_distance(intersection_points, intersection_distances, plot_context::PlotContext)
    trimmed_distances = Float64[]  # Initialize trimmed distances array

    # Iterate over intersection points to calculate trimmed lines
    for j in collect(2:2:lastindex(intersection_points))
        trimmed_line = [intersection_points[j-1], intersection_points[j]]
        trimmed_distance = distance_pointpoint(intersection_points[j-1], intersection_points[j])

        push!(trimmed_distances, trimmed_distance)  # Add trimmed distance

        # Plot if enabled
        if plot_context.plot
            plot_trimmed_line(plot_context.ax, trimmed_line, trimmed_distance)
        end
    end

    return sum(trimmed_distances)  # Calculate total distance
end
