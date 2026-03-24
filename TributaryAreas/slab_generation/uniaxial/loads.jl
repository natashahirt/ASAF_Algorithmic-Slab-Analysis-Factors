"""
    calculate_uniaxial_loads(cycle::Vector{Int64}, elements::Vector{<:Element}; vector_1d=[1.0, 0.0], spacing::Float64, w::Float64, fix_param=false, plot=false)

Calculate the point loads on each element in a 1D cycle.

# Arguments
- `cycle::Vector{Int64}`: Indices of nodes forming a cycle.
- `elements::Vector{<:Element}`: Elements to consider for load calculation.
- `vector_1d`: Direction vector for load application.
- `spacing::Float64`: Spacing between load application points.
- `w::Float64`: Load magnitude.
- `fix_param`: Flag to fix parameters during load calculation.
- `plot`: Flag to enable plotting.

# Returns
- `Vector{Asap.AbstractLoad}`: Calculated point loads for the elements.
"""
function calculate_uniaxial_loads(self::SlabAnalysisParams, cycle::Vector{Int64}, elements, w::Float64)
    startpoints, stiffnesses = get_start_and_end_coords(cycle, elements)
    @assert test_convexity(startpoints) "One of your cycles is not convex! (cycle nodes: $cycle)"

    cycle_edges = get_cycle_edges(startpoints)
    cycle_elements = get_cycle_elements(cycle, elements)
    point_loads = Asap.AbstractLoad[]

    for i in 1:lastindex(cycle_elements)
        edge = cycle_edges[i]
        element = cycle_elements[i]

        distances, params = calculate_distances_and_params_uniaxial(self, element, i, cycle_edges, cycle_elements, elements)

        @assert length(params) == length(distances) "Mismatch in parameters and distances."

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

function calculate_distances_and_params_uniaxial(self::SlabAnalysisParams, element, i, cycle_edges, cycle_elements, elements)

    params = Float64[]
    distances = Float64[]

    interp_x, interp_y, interp_params = interpolate_between_nodes(element.nodeStart, element.nodeEnd, self.perp ? self.perp_vector_1d : self.vector_1d, spacing=self.spacing)

    for j in 1:lastindex(interp_params)
        if interp_params[j] <= 0 || interp_params[j] >= 1
            continue
        end

        distance = find_opposite_edge(self, [interp_x[j], interp_y[j]], i, cycle_edges, cycle_elements, model_elements=elements)
        if isnothing(distance)
            continue
        else
            push!(distances, distance)
            push!(params, interp_params[j])
        end
    end

    return distances, params
end


"""
    get_slab_loads_uniaxial(self::SlabAnalysisParams; vector_1d=[1.0, 0.0], w::Float64, spacing::Float64, fix_param=false, plot=false, slab_sizer=:individual, perp::Bool=false)

Calculate and apply slab loads for a given slab analysis.

# Arguments
- `self::SlabAnalysisParams`: Parameters for slab analysis.
- `vector_1d`: Direction vector for load application.
- `w::Float64`: Load magnitude.
- `spacing::Float64`: Spacing between load application points.
- `fix_param`: Flag to fix parameters during load calculation.
- `plot`: Flag to enable plotting.
- `slab_sizer`: Method for slab sizing.
- `perp::Bool`: Flag for perpendicular plotting.

# Returns
- `SlabAnalysisParams`: Updated slab analysis parameters.
"""
function get_slab_loads_uniaxial(self::SlabAnalysisParams)

    slab_nodes = [node for node in self.model.nodes if node.id != :fixed]

    adjacency_dict, half_edges, slab_nodes = get_half_edges(slab_nodes, self.model.elements)
    cycles = get_cycles(adjacency_dict, half_edges)

    self.model.nodes = [slab_nodes; self.model.nodes[:fixed]]
    loads = Asap.AbstractLoad[]
    areas = Float64[]

    valid_cycles = filter_valid_cycles(self, cycles)

    max_spans, slab_depths = get_slab_depths(self, valid_cycles, self.model.elements[:beam])

    self.max_spans = max_spans
    self.slab_depths = slab_depths

    for (i, cycle) in enumerate(valid_cycles)
        point_loads = calculate_uniaxial_loads(self, cycle, self.model.elements[:beam], slab_depths[i])
        vertical_loads = [load.value[3] for load in point_loads]
        area = abs(sum(vertical_loads)) / slab_depths[i]
        append!(loads, point_loads)
        append!(areas, area)
    end

    area = sum(areas)
    self.areas = areas

    if self.plot_context.plot
        plot_node!(self.plot_context.ax, slab_nodes, color=:black, numbers=false)
        plot_elements!(self.plot_context.ax, self.model.elements[:beam], color=:black)
    end

    append!(self.model.loads, loads)
    Asap.solve!(self.model, reprocess=true)

    self.model, self.area = self.model, area
    return self
end