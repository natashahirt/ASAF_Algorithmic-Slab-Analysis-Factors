function identify_raster_points(slab_params::SlabAnalysisParams, resolution::Int=5)
    # Pre-allocate arrays for node positions
    node_positions_x = Vector{Float64}(undef, length(slab_params.model.nodes))
    node_positions_y = Vector{Float64}(undef, length(slab_params.model.nodes))
    
    # Batch process node positions
    for (i, node) in enumerate(slab_params.model.nodes)
        node_positions_x[i] = abs(node.position[1])
        node_positions_y[i] = abs(node.position[2])
    end
    
    # Calculate domain bounds
    max_x = maximum(node_positions_x)
    max_y = maximum(node_positions_y)

    # Use the smaller dimension to determine step size to ensure square cells
    step_size = min(max_x, max_y) / resolution
    
    # Calculate number of steps, ensuring whole numbers
    n_x = Int(floor(max_x / step_size))
    n_y = Int(floor(max_y / step_size))

    # Create raster points with uniform step size
    x_raster = range(0, max_x, length=n_x+1)
    y_raster = range(0, max_y, length=n_y+1)

    # Calculate centerpoints of cells
    centerpoints_x = [x + step_size/2 for x in x_raster[1:end-1]]
    centerpoints_y = [y + step_size/2 for y in y_raster[1:end-1]]

    # Pre-allocate centerpoint coordinates array
    n_points = n_x * n_y
    centerpoint_coords = Vector{Vector{Float64}}(undef, n_points)
    
    # Fill centerpoint coordinates more efficiently
    idx = 1
    for y in centerpoints_y, x in centerpoints_x
        centerpoint_coords[idx] = [x, y]
        idx += 1
    end

    return centerpoint_coords
end

function get_raster_df(slab_params::SlabAnalysisParams; resolution::Int=5)
    
    beams = slab_params.model.elements[:beam]

    # get the centerpoint coordinates
    centerpoint_coords = identify_raster_points(slab_params, resolution)
        
    # Pre-allocate vectors for DataFrame with estimated capacity
    estimated_size = length(beams) * length(centerpoint_coords) ÷ 4  # Rough estimate
    centerpoints = Vector{Vector{Float64}}(undef, estimated_size)
    centerpoint_ids = Vector{Int}(undef, estimated_size)
    beam_ids = Vector{Tuple}(undef, estimated_size)
    beam_idx = Vector{Int}(undef, estimated_size)
    distances = Vector{Float64}(undef, estimated_size)
    stiffnesses = Vector{Float64}(undef, estimated_size)
    t_values = Vector{Float64}(undef, estimated_size)
    x_values = Vector{Float64}(undef, estimated_size)
    areas = Vector{Float64}(undef, estimated_size)
    volumes = Vector{Float64}(undef, estimated_size)
    loads = Vector{Float64}(undef, estimated_size)
    
    idx = 1
    for (beam_i, beam) in enumerate(beams)
        beam_id = get_element_id(beam)
        beam_vector = beam.nodeEnd.position[1:2] - beam.nodeStart.position[1:2]
        beam_length = distance_pointpoint(beam.nodeStart, beam.nodeEnd)
        t_intervals = collect(0:slab_params.spacing:beam_length)/beam_length
        t_intervals[end] = 0.999
        t_intervals[1] = 0.001

        for (i, centerpoint) in enumerate(centerpoint_coords)
            point_vector = centerpoint - beam.nodeStart.position[1:2]
            t = project_to_line(beam_vector, point_vector)
            if 0 ≤ t ≤ 1
                distance = distance_pointline(beam.nodeStart, beam.nodeEnd, centerpoint)
                t_snapped = t_intervals[argmin(abs.(t_intervals .- t))]
                
                # Resize arrays if needed
                if idx > estimated_size
                    new_size = estimated_size * 2
                    resize!(centerpoints, new_size)
                    resize!(centerpoint_ids, new_size)
                    resize!(beam_ids, new_size)
                    resize!(beam_idx, new_size)
                    resize!(distances, new_size)
                    resize!(stiffnesses, new_size)
                    resize!(t_values, new_size)
                    resize!(x_values, new_size)
                    resize!(areas, new_size)
                    resize!(volumes, new_size)
                    resize!(loads, new_size)
                    resize!(beam_idx, new_size)
                    estimated_size = new_size
                end
                
                # Store values in pre-allocated arrays
                centerpoints[idx] = centerpoint
                centerpoint_ids[idx] = i
                beam_ids[idx] = beam_id
                beam_idx[idx] = beam_i
                distances[idx] = distance
                t_values[idx] = t_snapped
                x_values[idx] = t_snapped * beam_length
                idx += 1
            end
        end
    end
    
    # Trim arrays to actual size and create DataFrame
    actual_size = idx - 1
    slab_params.raster_df = DataFrame(
        centerpoint_coords=view(centerpoints, 1:actual_size),
        centerpoint_id=view(centerpoint_ids, 1:actual_size),
        beam_id=view(beam_ids, 1:actual_size),
        beam_idx=view(beam_idx, 1:actual_size),
        distance=view(distances, 1:actual_size),
        stiffnesses=zeros(actual_size),
        t=view(t_values, 1:actual_size),
        x=view(x_values, 1:actual_size),
        probability=zeros(actual_size),
        areas=zeros(actual_size),
        volumes=zeros(actual_size),
        loads=zeros(actual_size),
    )

    return slab_params

end

function get_load_probabilities(slab_params::SlabAnalysisParams)
    """
    Calculate load probabilities for each beam based on stiffness.
    """
    # Pre-compute beam lookup dictionary for faster access
    beam_lookup = Dict(id => beam for (id, beam) in pairs(slab_params.element_id_lookup_df))
    
    # Process each group of centerpoints
    for (i, group) in enumerate(groupby(slab_params.raster_df, :centerpoint_id))        
        # Vectorized stiffness calculation
        beams = [beam_lookup[id] for id in group.beam_id]
        group.stiffnesses .= get_equivalent_spring_constant.(beams, group.distance)
        
        # Vectorized probability calculation
        total_stiffness = sum(group.stiffnesses)
        group.probability .= group.stiffnesses ./ total_stiffness
    end

    return slab_params
end

function get_equivalent_spring_constant(beam::Asap.Element, distance::Float64)
    slab_stiffness = distance == 0 ? 1 : 1 / distance^3
    beam_stiffness = 1 / (beam.section.Ix / beam.length^3)
    return slab_stiffness + beam_stiffness
end

function get_raster_loads(slab_params::SlabAnalysisParams)
    """
    Calculate and apply tributary loads to the model based on raster analysis.
    """
    # Pre-compute coordinates and spacing
    x_coords = sort!(unique([row.centerpoint_coords[1] for row in eachrow(slab_params.raster_df)]))
    y_coords = sort!(unique([row.centerpoint_coords[2] for row in eachrow(slab_params.raster_df)]))
    
    tributary_area = minimum(diff(x_coords)) * minimum(diff(y_coords))
    slab_depth = slab_params.slab_depths[1]
    println("Slab depth: $(round(slab_depth, digits=2)) m")
    
    # Pre-allocate output arrays with estimated size
    n_beams = length(unique(slab_params.raster_df.beam_id))
    n_t_intervals = length(unique(slab_params.raster_df.t))
    estimated_size = n_beams * n_t_intervals
    
    loads = Vector{Asap.AbstractLoad}(undef, estimated_size)
    load_areas = Vector{Float64}(undef, estimated_size)
    load_volumes = Vector{Float64}(undef, estimated_size)

    # Process each beam group
    idx = 1
    beam_groups = groupby(slab_params.raster_df, :beam_id)
    for group in beam_groups
        element = slab_params.element_id_lookup_df[group.beam_id[1]]
        group.areas = group.probability .* tributary_area
        group.volumes = group.areas .* slab_depth

        # Process each t-interval
        for t_interval in unique(group.t)
            t_mask = group.t .== t_interval
            t_area = sum(view(group.areas, t_mask))
            t_volume = sum(view(group.volumes, t_mask))
            
            loads[idx] = Asap.PointLoad(element, t_interval, [0, 0, t_volume])
            load_areas[idx] = t_area
            load_volumes[idx] = t_volume
            idx += 1
        end
    end
    
    # Trim arrays to actual size
    actual_size = idx - 1
    resize!(loads, actual_size)
    resize!(load_areas, actual_size)
    resize!(load_volumes, actual_size)
    
    # Update slab_params
    slab_params.model.loads = loads
    slab_params.load_areas = load_areas
    slab_params.load_volumes = load_volumes
    
    n_centerpoints = length(unique(slab_params.raster_df.centerpoint_id))
    slab_params.areas = fill(tributary_area, n_centerpoints)
    
    println("Number of model loads: $(length(loads))")
    println("Number of load areas: $(length(load_areas))")
    
    return slab_params
end