"""
    get_start_and_end_coords(cycle::Vector{Int64}, elements::Vector{<:Element})

Given a cycle and a list of elements, this function returns the coordinates of the start nodes
and the stiffness ratios for each element in the cycle.
"""
function get_start_and_end_coords(cycle::Vector{Int64}, elements::Vector{<:Element})
    cycle_elements = get_cycle_elements(cycle, elements)
    cycle_node_coords = Vector[]
    cycle_stiffnesses = Float64[]

    # Iterate over each element in the cycle to determine start and end coordinates
    for i in 1:lastindex(cycle_elements)
        element1 = cycle_elements[mod1(i-1, lastindex(cycle_elements))]
        element2 = cycle_elements[i]

        # Determine the orientation of the nodes in the element
        if element1.nodeStart.nodeID in Asap.nodeids(element2)
            node_o = element1.nodeStart
            node_a = element1.nodeEnd
        else
            node_a = element1.nodeStart
            node_o = element1.nodeEnd
        end

        # Store the coordinates of the start node
        push!(cycle_node_coords, [node_o.position[1], node_o.position[2]])
        
        # Calculate and store the stiffness ratio
        stiffness_ratio = get_stiffness_ratio(element1, element2)
        push!(cycle_stiffnesses, 1 - stiffness_ratio)
    end

    return cycle_node_coords, cycle_stiffnesses
end

"""
    get_cycle_elements(cycle::Vector{Int64}, elements::Vector{<:Element})

Returns a vector of elements that form the given cycle.
"""
function get_cycle_elements(cycle::Vector{Int64}, elements::Union{Vector{Element{T}}, Vector{Element}}) where T <: Asap.Release
    cycle_elements = Element[]

    # Iterate over each node in the cycle to find corresponding elements
    for i in 1:length(cycle)
        element = get_element_from_edge([cycle[i], cycle[mod1(i+1, length(cycle))]], elements)
                
        # Add the found element to the cycle elements
        if !isnothing(element)
            push!(cycle_elements, element)
        else
            # Log the edge if no element is found
            println("no element found for edge: $([cycle[i], cycle[mod1(i+1, length(cycle))]])")
            return nothing
        end
    end

    return cycle_elements
end

"""
    get_element_from_edge(edge::Vector{Int64}, elements::Vector{<:Element})

Finds and returns the element that corresponds to the given edge.
"""
function get_element_from_edge(edge::Vector{Int64}, elements::Union{Vector{Element{T}}, Vector{Element}}) where T <: Asap.Release
    # Search for the element that contains both nodes of the edge
    for element in elements
        if edge[1] in Asap.nodeids(element) && edge[2] in Asap.nodeids(element)
            return element
        end
    end
    return nothing
end

"""
    get_all_edges(elements::Vector{<:Element})

Returns a vector of all edges in the given elements.
"""
function get_all_edges(elements::Vector{<:Element})
    edges = Vector[]

    # Collect start and end coordinates for each element
    for element in elements
        push!(edges, [[element.nodeStart.position[1], element.nodeStart.position[2]], 
                      [element.nodeEnd.position[1], element.nodeEnd.position[2]]])
    end

    return edges
end

"""
    get_cycle_edges(coordinates::Vector{Vector})

Returns a vector of edges that form a cycle from the given coordinates.
"""
function get_cycle_edges(coordinates::Vector{Vector})
    cycle_edges = Vector[]

    # Create edges by connecting consecutive coordinates
    for i in 1:length(coordinates)
        push!(cycle_edges, [coordinates[i], coordinates[mod1(i+1, length(coordinates))]])
    end

    return cycle_edges
end

"""
    interpolate_between_nodes(a::Node, b::Node, n::Int64)

Interpolates between two nodes `a` and `b` into `n` segments.
"""
function interpolate_between_nodes(a::Node, b::Node, n::Int64)
    a = [a.position[1], a.position[2]]
    b = [b.position[1], b.position[2]]
    
    return interpolate_between_nodes(a, b, n)
end

"""
    interpolate_between_nodes(element::Element, n::Int64)

Interpolates between the start and end nodes of an element into `n` segments.
"""
function interpolate_between_nodes(element::Element, n::Int64)
    nodeStart = element.nodeStart
    nodeEnd = element.nodeEnd

    return interpolate_between_nodes(nodeStart, nodeEnd, n)
end

"""
    interpolate_between_nodes(a::Node, b::Node, vector_1d::Vector{Float64}; spacing=0.1)

Interpolates between two nodes `a` and `b` using a given 1D vector and spacing.
"""
function interpolate_between_nodes(a::Node, b::Node, vector_1d::Vector{Float64}; spacing=0.1)
    vector_ab = [b.position[1] - a.position[1], b.position[2] - a.position[2]]
    distance_ab = sqrt(vector_ab[1]^2 + vector_ab[2]^2)

    # Determine interpolation method based on vector alignment
    if crossproduct(vector_1d, vector_ab) == 0 # perpendicular
        x_series, y_series, positions = [], [], [] # no interpolation needed
    elseif dotproduct(vector_1d, vector_ab) == 0 # parallel
        n = Int(floor(distance_ab / spacing))
        x_series, y_series, positions = interpolate_between_nodes(a, b, n)
    else # vectors are at an angle
        perp_vector = [vector_1d[2], -vector_1d[1]]
        scaled_vector = unit_vector(perp_vector) * spacing

        scalar_projection = (scaled_vector[1] * vector_ab[1] + scaled_vector[2] * vector_ab[2]) / distance_ab

        if scalar_projection < 0
            scaled_vector *= -1
            scalar_projection *= -1
        end

        scaled_vector_endpoint = [a.position[1], a.position[2]] + scaled_vector
        new_point, distance = ray_ray_intersect(scaled_vector_endpoint, vector_1d, [a.position[1], a.position[2]], vector_ab)
        t = distance_pointpoint([a.position[1], a.position[2]], new_point) / distance_ab
        t_current = t
        vector_projection = [new_point[1] - a.position[1], new_point[2] - a.position[2]]
        
        x_series = Float64[]
        y_series = Float64[]
        positions = Float64[]
        
        # Generate interpolated points along the vector
        while t_current < 1
            push!(x_series, new_point[1])
            push!(y_series, new_point[2])
            push!(positions, t_current)

            new_point = new_point + vector_projection
            t_current += t
        end
    end

    return x_series, y_series, positions
end

"""
    interpolate_between_nodes(a::Vector{Float64}, b::Vector{Float64}, n::Int64)

Interpolates between two points `a` and `b` into `n` segments.
"""
function interpolate_between_nodes(a::Vector{Float64}, b::Vector{Float64}, n::Int64)
    x1, x2, y1, y2 = a[1], b[1], a[2], b[2]
    x_spacing = (x2 - x1) / n
    y_spacing = (y2 - y1) / n
    param_spacing = 1 / n

    # Generate series of interpolated x, y coordinates and parameter positions
    x_series = [i * x_spacing + x1 - x_spacing / 2 for i in 1:n]
    y_series = [i * y_spacing + y1 - y_spacing / 2 for i in 1:n]
    positions = [i * param_spacing - param_spacing / 2 for i in 1:n]

    return x_series, y_series, positions
end

"""
    get_max_diagonal_span_2d(startpoints::Vector{Vector}; method::Symbol=:corners, vector_1d::Vector=[1.,0.])

Calculates the maximum diagonal span of a 2D shape defined by `startpoints` using the specified method.
"""
function get_max_diagonal_span_2d(self::SlabAnalysisParams, startpoints::Vector{Vector}; method::Symbol=:corners, vector_1d::Vector=[1.,0.])
    max_distance = 0.

    if method == :corners
        # Take the maximum distance between any set of corners. Most conservative.
        for i in 1:lastindex(startpoints)-1
            for j in i+1:lastindex(startpoints)
                distance = distance_pointpoint(startpoints[i], startpoints[j])
                if distance > max_distance
                    max_distance = distance
                end
            end
        end

    elseif method == :bisector
        # For each corner, find the bisector's intersection with the polygon and take that distance.
        spans = Float64[]
        edge_list = get_cycle_edges(startpoints)

        # Calculate bisector distances for each corner
        for i in 1:lastindex(startpoints)
            node_a = startpoints[mod1(i-1, lastindex(startpoints))]
            node_o = startpoints[i]
            node_b = startpoints[mod1(i+1, lastindex(startpoints))]

            edge_a = [node_a, node_o]
            edge_b = [node_o, node_b]

            if check_collinearity(edge_a, edge_b) == true
                edge_vector = node_b - node_a
                bisector_vector = [-edge_vector[2], edge_vector[1]]

                target_point, distance = find_opposite_edge_orthogonal(self, node_o, edge_list, vector_1d=bisector_vector, plot=false)

                if distance == 0
                    target_point, distance = find_opposite_edge_orthogonal(self, node_o, edge_list, vector_1d=-bisector_vector, plot=false)
                end
            else
                bisector = line_line_bisector(node_o, node_a, node_b, param=0.5)

                if isnan(bisector[1]) || isnan(bisector[2])
                    continue
                else
                    bisector_vector = bisector - node_o
                end

                target_point, distance = find_opposite_edge_orthogonal(self, node_o, edge_list, vector_1d=bisector_vector, plot=false)
            end

            if !isnan(distance) && distance != 0 && target_point != node_a && target_point != node_b
                push!(spans, distance)
            end
        end

        max_distance = maximum(spans)

    elseif method == :orthogonal
        # For each edge, find its normal vector and intersect with as many corners as possible. Take the maximum.
        spans = []
        edge_list = get_cycle_edges(startpoints)

        # Calculate orthogonal distances for each edge
        for edge in edge_list
            edge_vector = edge[2] - edge[1]

            for i in 1:lastindex(startpoints)
                node_a = startpoints[mod1(i-1, lastindex(startpoints))]
                node_o = startpoints[i]
                node_b = startpoints[mod1(i+1, lastindex(startpoints))]

                target_point, distance = find_opposite_edge_orthogonal(self, node_o, edge_list, vector_1d=self.perp_vector_1d, plot=false)

                if !isnan(distance) && distance != 0 && target_point != node_a && target_point != node_b
                    push!(spans, distance)
                end
            end
        end

        if isempty(spans)
            spans = [distance_pointpoint(edge[1], edge[2]) for edge in edge_list]
        end

        max_distance = maximum(spans)

    elseif method == :orthonormal
        # Using vector 1d, take the maximum span.
        perp_vector_1d = [vector_1d[2], -vector_1d[1]]

        cycle_edges = get_cycle_edges(startpoints)

        # Calculate maximum spans in both parallel and perpendicular directions
        max_span_parallel = get_max_diagonal_span_1d(cycle_edges, vector_1d)
        max_span_perp = get_max_diagonal_span_1d(cycle_edges, perp_vector_1d)

        max_distance = maximum([max_span_parallel, max_span_perp])
    end

    return max_distance
end

"""
    get_max_diagonal_span_1d(edge_list::Vector{Vector}, vector_1d::Vector{Float64})

Calculates the maximum diagonal span of a 1D shape defined by `edge_list` using the given vector.
"""
function get_max_diagonal_span_1d(edge_list::Vector{Vector}, vector_1d::Vector{Float64})
    max_distance = 0.

    # Iterate over each edge to find the maximum span
    for i in 1:lastindex(edge_list)
        point = edge_list[i][1] # edge startpoint
        
        for j in 1:lastindex(edge_list)
            target_point, distance = ray_line_intersect(point, vector_1d, edge_list[j], param=1.)

            if target_point == point
                continue
            end

            if distance > max_distance
                max_distance = distance
            end
        end
    end

    return max_distance
end

"""
    get_indices_nonzero_edges(edges; tol::Float64=0.01)

Returns the indices of edges that exceed a given tolerance length.
"""
function get_indices_nonzero_edges(edges; tol::Float64=0.01)
    i_nonzero_edges = Int64[]

    # Identify edges with length greater than the tolerance
    for (i, edge) in enumerate(edges)
        if get_length(edge) > tol
            push!(i_nonzero_edges, i)
        end
    end

    return i_nonzero_edges
end

"""
    scale_slab(self::SlabAnalysisParams, scale=1.; dimensions=2)

Scales the slab model by a given factor. This function adjusts the positions of nodes in the slab model and the spacing parameter based on the specified scale. It is a destructive operation that modifies the original model.

# Arguments
- `self::SlabAnalysisParams`: The slab analysis parameters containing the model to be scaled.
- `scale=1.`: The scaling factor. The function uses the square root of this value to normalize the area.
- `dimensions=2`: The number of dimensions to scale. Defaults to 2, but can be set to 3 for 3D scaling.

# Returns
- `self`: The modified `SlabAnalysisParams` object with updated node positions and spacing.
"""
function scale_slab(self::SlabAnalysisParams, scale=1.; dimensions=2)
    # Calculate the scaling factor for length from the area scaling factor
    scale = sqrt(scale)

    # Scale the position of each node in the model
    for node in self.model.nodes
        node.position[1] *= scale
        node.position[2] *= scale

        # Scale the third dimension if specified
        if dimensions == 3
            node.position[3] *= scale
        end
    end

    # Scale the spacing parameter
    self.spacing *= scale

    return self
end