function plot_centerpoint_lines(slab_params::SlabAnalysisParams; centerpoint_id::Int=0, reset_plot::Bool=true, text::Bool=true, ax::Union{Nothing, Axis}=nothing)
    
    if (isnothing(slab_params.plot_context.fig) || reset_plot) && isnothing(ax)
        slab_params.plot_context = setup_plot(slab_params.plot_context.plot; size=(1200, 800))
        plot_model(slab_params.model, plot_context=slab_params.plot_context)
    end

    if centerpoint_id > 0
        # Filter data for the specific centerpoint_id
        df = filter(row -> row.centerpoint_id == centerpoint_id, slab_params.raster_df)
    else
        df = slab_params.raster_df
    end
    
    # Pre-allocate arrays for better performance
    n_points = nrow(df)
    centerpoint_xs = Vector{Float64}(undef, n_points)
    centerpoint_ys = Vector{Float64}(undef, n_points)
    line_xs = Vector{Vector{Float64}}(undef, n_points)
    line_ys = Vector{Vector{Float64}}(undef, n_points)
    distances = Vector{Float64}(undef, n_points)
    
    # Batch process the data
    for (i, row) in enumerate(eachrow(df))
        beam_id = row.beam_id
        beam = slab_params.element_id_lookup_df[beam_id]
        nodestart = beam.nodeStart
        nodeend = beam.nodeEnd
        beam_vector = nodeend.position[1:2] - nodestart.position[1:2]
        beam_intersect_point = nodestart.position[1:2] + row.t * beam_vector
        
        # Store data in pre-allocated arrays
        centerpoint_xs[i] = row.centerpoint_coords[1]
        centerpoint_ys[i] = row.centerpoint_coords[2]
        line_xs[i] = [row.centerpoint_coords[1], beam_intersect_point[1]]
        line_ys[i] = [row.centerpoint_coords[2], beam_intersect_point[2]]
        distances[i] = row.distance
    end
        
    # Batch plot operations
    scatter!(slab_params.plot_context.ax, centerpoint_xs, centerpoint_ys, color="red")
    
    # Plot lines only for the selected centerpoint
    for i in 1:n_points
        if df.probability[i] > 0.01
            lines!(slab_params.plot_context.ax, line_xs[i], line_ys[i], 
                color=slab_params.raster_df.probability[i], colormap=Reverse(:acton), 
                colorrange=(0, 1))
            if text
                text!(slab_params.plot_context.ax, 
                    mean(line_xs[i]), mean(line_ys[i]), 
                    text=string(round(df.probability[i], digits=2)),
                    align=(:center, :center),
                    color=:black,
                    fontsize=12)
            end
        end
    end
    
    return slab_params.plot_context.fig
end

"""
Create an animation showing probabilities for each centerpoint sequentially.

Parameters
----------
slab_params : SlabAnalysisParams
    The slab analysis parameters containing the model and probability data
kwargs : Dict
    Optional keyword arguments:
    - framerate : Int = 2
    - output_file : String = "centerpoint_probabilities.gif"

Returns
-------
String
    Path to the generated animation file
"""
function animate_centerpoint_probabilities(slab_params::SlabAnalysisParams; 
    framerate::Int=2, 
    output_file::String="centerpoint_probabilities_925_clean.gif")
    
    unique_centerpoints = sort(unique(slab_params.raster_df.centerpoint_id))
    
    # Create figure and axis once
    fig = Figure(size=(1200, 800))
    ax = Axis(fig[1,1])
    
    # Create observable for current centerpoint
    current_point = Observable(1)
    
    # Plot base model once
    plot_model(slab_params.model, plot_context=slab_params.plot_context, ax=ax)
    
    # Create observables for the visualization
    df = slab_params.raster_df
    
    # Create the visualization that updates based on current_point
    @lift begin
        # Filter data for the specific centerpoint_id
        point_df = filter(row -> row.centerpoint_id == unique_centerpoints[$current_point], df)
        
        # Plot centerpoints and lines
        for row in eachrow(point_df)
            beam = slab_params.element_id_lookup_df[row.beam_id]
            beam_vector = beam.nodeEnd.position[1:2] - beam.nodeStart.position[1:2]
            beam_intersect_point = beam.nodeStart.position[1:2] + row.t * beam_vector
            
            # Plot centerpoint
            scatter!(ax, [row.centerpoint_coords[1]], [row.centerpoint_coords[2]], color="red")
            
            # Add probability text if needed
            if row.probability > 0.01
                # Plot connection line
                lines!(ax, 
                    [row.centerpoint_coords[1], beam_intersect_point[1]], 
                    [row.centerpoint_coords[2], beam_intersect_point[2]],
                    color=row.probability,
                    colormap=Reverse(:acton),
                    colorrange=(0, 1))
                """text!(ax,
                    mean([row.centerpoint_coords[1], beam_intersect_point[1]]),
                    mean([row.centerpoint_coords[2], beam_intersect_point[2]]),
                    text=string(round(row.probability, digits=2)),
                    align=(:center, :center),
                    color=:black,
                    fontsize=12)"""
            end
        end
    end
    
    # Create animation
    println("Generating animation for $(length(unique_centerpoints)) frames...")
    record(fig, output_file, 1:length(unique_centerpoints); framerate=framerate) do frame
        current_point[] = frame
        
        # Progress reporting
        if frame % 5 == 0
            println("Processing frame $frame/$(length(unique_centerpoints))")
        end
    end
    
    return output_file
end
