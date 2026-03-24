"""
    generate_trapezium(h::Real, b_lower::Real, b_upper::Real; span=false, i_columns = [1,2,3,4], plot=false)

Generates a trapezoidal slab geometry with a lower base `b_lower`, upper base `b_upper`, and height `h`.

# Arguments
- `h::Real`: Height of the trapezium.
- `b_lower::Real`: Lower base length.
- `b_upper::Real`: Upper base length.
- `span::Bool`: If true, the trapezium spans from one side to the other.
- `i_columns`: Indices of columns to be included in the model.
- `plot::Bool`: Whether to plot the model. Defaults to false.

# Returns
- `model`: A preprocessed model with nodes and elements.
"""
function generate_trapezium(h::Real, b_lower::Real, b_upper::Real; span=false, i_columns = [1,2,3,4], plot=false)
    # Calculate the offset for the upper base
    offset = (b_lower - b_upper) / 2

    if span == false
        # Define positions for non-spanning trapezium
        positions = [
            [0., 0., 0.],
            [offset, h, 0.],
            [b_upper + offset, h, 0.],
            [b_lower, 0., 0.]
        ]

        # Define connectivity for non-spanning trapezium
        connectivity = [[1, 2], [2, 3], [3, 4], [4, 1]]

    else
        # Define positions for spanning trapezium
        positions = [
            [0., 0., 0.],
            [offset, h, 0.],
            [b_upper + offset, h, 0.],
            [b_lower, 0., 0.],
            [b_upper + offset, 0., 0.],
            [offset, 0., 0.]
        ]

        # Define connectivity for spanning trapezium
        connectivity = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1], [2, 6], [3, 5]]
    end

    # Create nodes and elements
    nodes = Node[Node(position, :free) for position in positions]
    elements = connect_nodes(nodes, connectivity)

    # Assign IDs to nodes and elements
    for node in nodes node.id = :slab end
    for element in elements element.id = :beam end

    # Create and preprocess the model
    model = Asap.Model(nodes, elements, Asap.AbstractLoad[])
    model = preprocess_model(model, i_columns=i_columns)

    # Plot the model if requested
    if plot
        plot_model(model)
    end

    return model
end
