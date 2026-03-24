"""
    generate_from_json(dict::Dict; plot=true, sections::Vector=[], drawn::Bool=false)

Generate a structural model from a JSON-like dictionary.

# Arguments
- `dict::Dict`: A dictionary containing the structure's data.
- `sections::Vector`: A vector of sections to be used in the model. Default is an empty vector.
- `drawn::Bool`: Indicates if the elements have been drawn. Default is `false`.
- `plot::Bool`: Whether to plot the model. Default is `true`.

# Returns
- `model`: An `Asap.Model` object representing the structural model.
"""
function generate_from_json(dict::Dict; sections::Vector{Section}=Section[], drawn::Bool=false, plot=false) 
    # Extract positions and connectivities from the dictionary
    positions_dict = dict["positions"]
    connectivities_dict = dict["connectivities"]

    # Convert positions and connectivities to arrays
    # Convert positions to array and find minimum x,y,z coordinates
    positions = [[position["x"], position["y"], position["z"]] for position in positions_dict]
    min_x = minimum(p[1] for p in positions)
    min_y = minimum(p[2] for p in positions) 
    min_z = minimum(p[3] for p in positions)
    
    # Subtract minimum coordinates to start at origin
    positions = [[p[1] - min_x, p[2] - min_y, p[3] - min_z] for p in positions]
    connectivities = [([Int64(connect["to"]), Int64(connect["from"])]) for connect in connectivities_dict]
    i_columns = [Int64(i) for i in dict["columns_i"]]
    # Check if wall_i exists in dictionary, if not use empty array
    i_wall = haskey(dict, "wall_i") ? [Int64(i) for i in dict["wall_i"]] : Int64[]
    i_holes = haskey(dict, "holes_i") ? vcat([Int64[j for j in i["hole"]] for i in dict["holes_i"]]...) : Int64[]
    i_perimeter = haskey(dict, "perimeter_i") ? [Int64(i) for i in dict["perimeter_i"]] : Int64[]
    
    nodes = Node[Node(position, :free) for position in positions]

    column_set = Set(i_columns)
    wall_set = Set(i_wall)
    for (i, node) in enumerate(nodes)
        if i in column_set
            node.id = :column
            node.dof = [true, true, true, false, false, false]
        elseif i in wall_set
            node.id = :wall
            node.dof = [false, false, false, false, false, false]
        else
            node.id = :slab
            node.dof = [true, true, true, false, false, false]
        end
    end

    # Determine element IDs if available and drawn is true
    element_ids = haskey(connectivities_dict[1], "id") && drawn == true ? 
        [Int64(connect["id"]) for connect in connectivities_dict] : Int64[]

    material = material_dict[:in]
    if isempty(sections)
        default_section = toASAPframe_W(W_imperial(W_dict_imperial["W8X35"]), material.E, material.G; convert=true)::Asap.Section # use mm so that there's no factoring
        sections = [default_section for _ in 1:lastindex(connectivities)]
    end

    default_column_section = toASAPframe_W(W_imperial(W_dict_imperial["W43X335"]), material.E, material.G; convert=true)::Asap.Section # use mm so that there's no factoring

    releases, connectivities = get_beam_releases(nodes, connectivities, element_ids=element_ids, i_holes=i_holes)
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

    if haskey(dict, "holes_i")
        i_holes = [Int64[j for j in i["hole"]] for i in dict["holes_i"]]
    else
        i_holes = Int64[]
    end

    # Make geometry data Dict
    type_information = Dict(
        "i_perimeter" => i_perimeter,
        "i_holes" => i_holes,
        "element_ids" => element_ids
    )

    # Plot the model if requested
    if plot
        plot_context = setup_plot()
        plot_model(model, plot_context=plot_context, numbers=true, type_information=type_information)
        display(plot_context.fig)
    end
    
    return model, type_information
end
