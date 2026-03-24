"""
    get_element_loads_isotropic(startpoints, cycle_elements; w, spacing, fix_param, plot=false, orthogonal, vector_1d=[1.,0.], model_elements)

Calculate the point loads on 2D elements based on the given parameters.

# Arguments
- `startpoints::Vector{Vector}`: Starting points of the elements.
- `cycle_elements::Vector{<:Element}`: Elements in the cycle.
- `w::Real`: Load weight.
- `spacing::Real`: Spacing between loads.
- `fix_param::Bool`: Whether to fix parameters.
- `plot::Bool`: Whether to plot the results.
- `orthogonal::Bool`: Whether the load is orthogonal.
- `vector_1d::Vector`: Direction vector for 1D elements.
- `model_elements::Vector{<:Element}`: Model elements.

# Returns
- `point_loads`: A vector of point loads.
"""
function get_element_loads_isotropic(self::SlabAnalysisParams, startpoints::Vector{Vector}, cycle_elements::Vector{<:Element}, w::Real, tol::Real=0.1)
    
    interior_edges = get_cycle_edges(startpoints)
    edge_levels = [interior_edges]
    i_nonzero_edges = get_indices_nonzero_edges(interior_edges, tol=tol)
    
    i = 1

    while length(i_nonzero_edges) > 2 && i < 100 # Ensure at least a triangle
        bisectors, intersections, interior_edges = get_interior_polygon(interior_edges, i_nonzero_edges, self.fix_param)
        push!(edge_levels, interior_edges)

        i_nonzero_edges = get_indices_nonzero_edges(interior_edges, tol=tol)

        if self.record_tributaries
            # Convert edges to point series
            interior_polygon = Vector[]
            if !isempty(i_nonzero_edges)
                # Add first point
                
                first_edge = interior_edges[i_nonzero_edges[1]]
                push!(interior_polygon, first_edge[1])
            
                # Add remaining points in order
                for idx in i_nonzero_edges
                    if !isapprox(interior_edges[idx][2], first_edge[1], atol=tol) && 
                       !any(p -> isapprox(p, interior_edges[idx][2], atol=tol), interior_polygon)
                        push!(interior_polygon, interior_edges[idx][2])
                    end
                end

                if test_convexity(interior_polygon)
                    push!(self.trib_dictionary["interior_polygons"], interior_polygon)
                end
            end
        end

        i += 1
    end

    sorted_by_edge = reorganize_edge_levels(edge_levels)
    nodes_list = convert_edges_to_nodes(sorted_by_edge)

    if self.record_tributaries
        for (i, edge_nodes) in enumerate(nodes_list)
            element_id = (cycle_elements[i].nodeStart.nodeID, cycle_elements[i].nodeEnd.nodeID)
            # Remove duplicate nodes from edge_nodes
            unique_nodes = Vector{Vector{Float64}}()
            for node in edge_nodes
                if !any(n -> isapprox(n, node, atol=1e-10), unique_nodes)
                    push!(unique_nodes, node)
                end
            end
            if length(unique_nodes) == 3
                push!(self.trib_dictionary["interior_polygons"], [unique_nodes[2]])
            end
            push!(self.trib_dictionary[element_id]["trib_area"], unique_nodes)
        end
    end
    
    point_loads = calculate_isotropic_loads(self, cycle_elements, nodes_list, w)

    return point_loads
end

function reorganize_edge_levels(edge_levels)
    sorted_by_edge = [[] for _ in edge_levels[1]]
    for level_i in 1:lastindex(edge_levels)
        for j in 1:lastindex(edge_levels[level_i])
            push!(sorted_by_edge[j], edge_levels[level_i][j])
        end
    end
    return sorted_by_edge
end

function convert_edges_to_nodes(sorted_by_edge)
    nodes_list = []
    for n in 1:lastindex(sorted_by_edge)
        edge_list = sorted_by_edge[n]
        node_list = []

        for edge in edge_list
            push!(node_list, edge[1])
        end

        for edge in reverse(edge_list)
            push!(node_list, edge[2])
        end

        smooth_node_list = smooth_nodes(node_list)
        push!(nodes_list, smooth_node_list)
    end
    return nodes_list
end

function smooth_nodes(node_list)
    smooth_node_list = [node_list[1], node_list[2]]
    i = 1
    while i <= lastindex(node_list)
        if i == lastindex(node_list)
            push!(smooth_node_list, node_list[end])
        end

        if check_collinearity(smooth_node_list[end-1], smooth_node_list[end], node_list[i], tol=0.1)
            smooth_node_list[end] = node_list[i]
        else
            push!(smooth_node_list, node_list[i])
        end
        i += 1
    end
    return smooth_node_list
end

function calculate_isotropic_loads(self::SlabAnalysisParams, cycle_elements, nodes_list, w)
    point_loads = Asap.AbstractLoad[]

    for (i,element) in enumerate(cycle_elements)

        """node_list = nodes_list[i]

        edges = [[nodes_list[i][j-1], nodes_list[i][j]] for j in 2:lastindex(nodes_list[i]) if nodes_list[i][j-1] != nodes_list[i][j]]
        for edge in edges
            plot_line!(self.plot_context.ax, edge, color=:lightgrey)
        end"""

        if !self.record_tributaries
            if self.fix_param
                params, distances, vector = get_element_params(self, element, nodes_list[i])
            else
                params, distances = get_element_params_nonconvex(self,element, nodes_list[i])
            end
        else
            vector = nothing
            params, distances, vector = get_element_params(self, element, nodes_list[i])
            # Create or append to load dictionary for this element
            if !isnothing(vector)
                element_id = (element.nodeStart.nodeID, element.nodeEnd.nodeID)
                append!(self.trib_dictionary[element_id]["trib_lines"], [(d, p, vector) for (d, p) in zip(distances, params)])
            end
        end

        # Calculate areas and volumes for each strip
        load_areas = distances .* self.spacing
        load_volumes = load_areas .* w
        append!(self.load_areas, load_areas)
        append!(self.load_volumes, load_volumes)
        append!(self.load_widths, distances)
        
        # Create load vectors with vertical (z) component only
        load_vectors = [[0.0, 0.0, volume] for volume in load_volumes]
        append!(point_loads, [PointLoad(element, params[i], load_vectors[i]) for i in 1:lastindex(load_vectors)])
    end
    return point_loads
end

"""
    get_slab_loads_isotropic(self; w, spacing, fix_param=false, plot=false, slab_sizer=:individual, orthogonal=false, vector_1d=[1.,0.])

Calculate the slab loads for a 2D slab analysis.

# Arguments
- `self::SlabAnalysisParams`: The slab analysis parameters.
- `w::Real`: Load weight.
- `spacing::Real`: Spacing between loads.
- `fix_param::Bool`: Whether to fix parameters.
- `plot::Bool`: Whether to plot the results.
- `slab_sizer`: Slab sizing method.
- `orthogonal::Bool`: Whether the load is orthogonal.
- `vector_1d::Vector`: Direction vector for 1D elements.

# Returns
- `self`: Updated slab analysis parameters.
"""
function get_slab_loads_isotropic(self::SlabAnalysisParams)
    
    slab_nodes = [node for node in self.model.nodes if node.id != :fixed]

    adjacency_dict, half_edges, slab_nodes = get_half_edges(slab_nodes, self.model.elements)
    cycles = get_cycles(adjacency_dict, half_edges)
    self.record_tributaries ? self.trib_dictionary["cycles"] = cycles : nothing
    
    self = calculate_slab_loads(self, cycles)

    self.area = sum(self.areas)

    println("Number of model loads: $(length(self.model.loads))")

    Asap.solve!(self.model, reprocess=true)

    if self.plot_context.plot
        plot_model(self.model, plot_context=self.plot_context)
    end

    return self
end

function calculate_slab_loads(self::SlabAnalysisParams, cycles)
    
    self.model.loads = Asap.AbstractLoad[]
    self.areas = Float64[]
    valid_cycles = filter_valid_cycles(self, cycles)

    self.max_spans, self.slab_depths = get_slab_depths(self, valid_cycles, self.model.elements[:beam])

    if self.load_type == :determinate
    for (i, cycle) in enumerate(valid_cycles)
        startpoints, stiffnesses = get_start_and_end_coords(cycle, self.model.elements[:beam])
        cycle_elements = get_cycle_elements(cycle, self.model.elements[:beam])  
            
            point_loads = get_element_loads_isotropic(self, startpoints, cycle_elements, self.slab_depths[i])
            vertical_loads = [load.value[3] for load in point_loads]

            if isempty(vertical_loads)
                continue
            end

            area = abs(sum(vertical_loads)) / self.slab_depths[i]
            append!(self.model.loads, point_loads)
            append!(self.areas, area)
        end
    elseif self.load_type == :indeterminate
        self = calculate_slab_loads_indeterminate(self)
    end

    return self
end