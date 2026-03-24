"""
    generate_funky_test(i_columns=Int64[])

Generates a model with a predefined funky geometry.

# Arguments
- `i_columns::Vector{Int64}`: Indices of columns to be included in the model. Defaults to an empty vector.
- `plot::Bool`: Whether to plot the model. Defaults to false.

# Returns
- `model::Asap.Model`: The generated model with nodes and elements.
"""
function generate_funky_test(i_columns=Int64[]; plot=false)

    # Define node positions in 3D space
    positions = [
        [0., -1., 0], [0., 10, 0], [12., 10., 0], [5., 0, 0],
        [1.5, 2.5, 0], [1.5, 7.5, 0], [3.5, 7.5, 0], [3.5, 2.5, 0],
        [2., 5., 0], [-2., 2.5, 0], [-2., 7.2, 0], [2.5, -2., 0]
    ]

    # Define connectivity between nodes
    connectivities = [
        [2,11], [11,10], [10,1], [1,4], [4,3], [3,2], [2,6], [1,5],
        [8,4], [8,3], [3,7], [6,7], [7,8], [8,5], [6,5], [6,9],
        [7,9], [9,8], [5,9], [12,1], [12,4]
    ]

    # Create nodes and elements
    nodes = [Node(position, :free) for position in positions]
    elements = connect_nodes(nodes, connectivities)

    # Assign IDs to nodes and elements
    for node in nodes
        node.id = :slab
    end
    for element in elements
        element.id = :beam
    end

    # Create and preprocess the model
    model = Asap.Model(nodes, elements, Asap.AbstractLoad[])
    model = preprocess_model(model, i_columns=i_columns)

    # Plot the model if requested
    if plot
        plot_model(model)
    end

    return model

end
