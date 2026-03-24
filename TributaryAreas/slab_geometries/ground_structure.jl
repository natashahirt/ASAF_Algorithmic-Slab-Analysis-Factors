"""
    generate_nested_rectangle(h_outer, b_outer, h_inner, b_inner; rotation=0., i_columns=[1,2,3,4])

Generates a nested rectangle model with optional rotation.

# Arguments
- `h_outer::Real`: Height of the outer rectangle.
- `b_outer::Real`: Base of the outer rectangle.
- `h_inner::Real`: Height of the inner rectangle.
- `b_inner::Real`: Base of the inner rectangle.
- `rotation::Real`: Rotation angle in radians. Defaults to 0.
- `i_columns::Vector{Int64}`: Indices of columns to be included in the model. Defaults to [1,2,3,4].
- `plot::Bool`: Whether to plot the model. Defaults to false.

# Returns
- `model::Asap.Model`: The generated model with nodes and elements.
"""
function generate_ground_structure(h::Real, w::Real, resolution::Real; no_diagonal_connections::Bool=false, moment_release::Bool=false, plot::Bool=false, sections::Vector{Asap.Section}=Asap.Section[], i_wall=[], element_ids::Vector{Int64}=Int64[])

    max_x = w
    max_y = h

    # Use the smaller dimension to determine step size to ensure square cells
    step_size = min(max_x, max_y) / resolution
    
    # Calculate number of steps, ensuring whole numbers
    n_x = Int(floor(max_x / step_size))
    n_y = Int(floor(max_y / step_size))

    # Create raster points with uniform step size
    # Create grid points including edges
    x_raster = range(0, max_x, length=n_x+1)
    y_raster = range(0, max_y, length=n_y+1)

    # Create nodes and assign them as free
    positions = [[x, y, 0.0] for x in x_raster for y in y_raster]

    # Add corner nodes and get their indices
    corner_positions = [
        [0.0, 0.0, 0.0],     # Bottom left
        [0.0, max_y, 0.0],   # Top left  
        [max_x, max_y, 0.0], # Top right
        [max_x, 0.0, 0.0]    # Bottom right
    ]

    # Add corner nodes only if they are not already in positions
    for corner in corner_positions
        if !any(p -> p == corner, positions)
            push!(positions, corner)
        end
    end
    
    # Get indices of corner nodes
    i_columns = [findfirst(p -> p == corner, positions) for corner in corner_positions]

    nodes = Node[Node(position, :free) for position in positions]

    # Define connectivity for each node to its eight neighbors
    connectivities_set = Set{Tuple{Int64, Int64}}()
    for i in 1:length(x_raster)
        for j in 1:length(y_raster)
            current_index = (i - 1) * length(y_raster) + j
            if no_diagonal_connections
                neighbors = [
                    (i-1, j), (i, j+1),
                    (i+1, j), (i, j-1)
                ]
            else
                neighbors = [
                (i-1, j-1), (i-1, j), (i-1, j+1),
                (i, j-1),           (i, j+1),
                    (i+1, j-1), (i+1, j), (i+1, j+1)
                ]
            end
            for (ni, nj) in neighbors
                if ni >= 1 && ni <= length(x_raster) && nj >= 1 && nj <= length(y_raster)
                    neighbor_index = (ni - 1) * length(y_raster) + nj
                    # Ensure the smaller index is first to avoid duplicates
                    connection = (min(current_index, neighbor_index), max(current_index, neighbor_index))
                    push!(connectivities_set, connection)
                end
            end
        end
    end
    
    # Convert the set back to a vector
    connectivities = [[c[1], c[2]] for c in collect(connectivities_set)]

    # Create nodes and assign them to the slab
    for (i, node) in enumerate(nodes)
        if i in i_columns
            node.id = :column
            node.dof = [true, true, true, false, false, false] # welded rigid connections
        elseif i in i_wall
            node.id = :wall
            node.dof = [false, false, false, false, false, false] # no degrees of freedom (concrete)
        else
            node.id = :slab
            node.dof = [true, true, true, false, false, false] # welded rigid connections
        end
    end

    # Determine element IDs if available and drawn is true
    """element_ids = haskey(connectivities_dict[1], "id") && drawn == true ? 
        [Int64(connect["id"]) for connect in connectivities_dict] : Int64[]"""

    material = material_dict[:in]
    if isempty(sections)
        default_section = toASAPframe_W(W_imperial(W_dict_imperial["W8X35"]), material.E, material.G; convert=true)::Asap.Section # use mm so that there's no factoring
        sections = [default_section for _ in 1:lastindex(connectivities)]
    end

    default_column_section = toASAPframe_W(W_imperial(W_dict_imperial["W43X335"]), material.E, material.G; convert=true)::Asap.Section # use mm so that there's no factoring

    releases, connectivities = get_beam_releases(nodes, connectivities, element_ids=element_ids)
    
    if moment_release
        """for (i,connection) in enumerate(connectivities)
            node1 = nodes[connection[1]]
            node2 = nodes[connection[2]]
            node1_free_x = false
            node2_free_x = false
            node1_free_y = false
            node2_free_y = false
            if node1.position[1] == 0.0 || node1.position[1] == max_x
                node1_free_x = true
            end
            if node2.position[1] == 0.0 || node2.position[1] == max_x
                node2_free_x = true
            end
            if node1.position[2] == 0.0 || node1.position[2] == max_y
                node1_free_y = true
            end
            if node2.position[2] == 0.0 || node2.position[2] == max_y
                node2_free_y = true
            end
            if node1_free_x && node2_free_x || node1_free_y && node2_free_y
                releases[i] = :fixedfixed
            elseif node1_free_x && node2_free_y || node1_free_y && node2_free_x
                releases[i] = :freefree
            elseif node1_free_x || node1_free_y
                releases[i] = :freefixed
            elseif node2_free_x || node2_free_y
                releases[i] = :fixedfree
            end
        end"""

        """for (i,connection) in enumerate(connectivities)
            node1 = nodes[connection[1]]
            node2 = nodes[connection[2]]
            if abs(node1.position[1] - node2.position[1]) < 1e-3 && !(node1.position[1] == 0.0 || node1.position[1] == max_x)
                releases[i] = :freefree
            end
        end"""

        x_unique = sort(unique([node.position[1] for node in nodes]))
        y_unique = sort(unique([node.position[2] for node in nodes]))
        """for (i, connection) in enumerate(connectivities)
            node1 = nodes[connection[1]]
            node2 = nodes[connection[2]]
            if node1.position[1] in x_unique && (node2.position[1] == x_unique[mod1(findfirst(==(node1.position[1]), x_unique) + 1, length(x_unique))] || node2.position[1] == x_unique[mod1(findfirst(==(node1.position[1]), x_unique) - 1, length(x_unique))]) && node1.position[2] != node2.position[2]
                releases[i] = :freefree
            end
        end"""

        for (i, connection) in enumerate(connectivities)
            node1 = nodes[connection[1]]
            node2 = nodes[connection[2]]
            # Check if element goes diagonally left (x decreases and y changes)
            if node2.position[1] == x_unique[mod1(findfirst(==(node1.position[1]), x_unique) + 1, length(x_unique))] && 
               node2.position[2] == y_unique[mod1(findfirst(==(node1.position[2]), y_unique) + 1, length(y_unique))]
                releases[i] = :freefree
            end
        end
    end

    beam_elements = connect_nodes(nodes, connectivities, releases=releases, sections=sections)
    nodes, column_elements = connect_columns(nodes, i_columns=i_columns, section=default_column_section)

    # Assign IDs to elements
    for element in beam_elements
        element.id = :beam
    end
    for element in column_elements
        element.id = :column
    end

    # Print the number of nodes in each category
    println("Number of slab nodes: $(length(nodes[:slab]))")
    println("Number of column nodes: $(length(nodes[:column]))")
    println("Number of wall nodes: $(length(nodes[:wall]))")
    println("Number of fixed nodes: $(length(nodes[:fixed]))")

    # Combine elements and create the model
    elements = [beam_elements; column_elements]
    model = Asap.Model(nodes, elements, Asap.AbstractLoad[]);

    # Plot the model if requested
    if plot
        plot_context = setup_plot()
        plot_model(model, plot_context=plot_context, numbers=true)
        display(plot_context.fig)
    end

    return model
end
