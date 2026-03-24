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
function generate_nested_rectangle(h_outer::Real, b_outer::Real, h_inner::Real, b_inner::Real; rotation=0., i_columns=[1,2,3,4], plot=false)

    # Define initial positions of the outer and inner rectangles
    positions = [
        [0., 0., 0.], [0., h_outer, 0.], [b_outer, h_outer, 0.], [b_outer, 0., 0.],
        [(b_outer-b_inner)/2, (h_outer-h_inner)/2., 0.],
        [(b_outer-b_inner)/2, (h_outer-h_inner)/2. + h_inner, 0.],
        [(b_outer-b_inner)/2 + b_inner, (h_outer-h_inner)/2. + h_inner, 0.],
        [(b_outer-b_inner)/2 + b_inner, (h_outer-h_inner)/2., 0.]
    ]

    # Define connectivity between nodes
    connectivity = [[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5],[1,5],[2,6],[3,7],[4,8]]

    # Rotate inner rectangle if rotation is specified
    if rotation != 0
        centroid = [b_outer / 2, h_outer / 2]  # Calculate centroid of the outer rectangle
        for i in 5:8
            positions[i] = rotate_point(centroid, positions[i], rotation)  # Rotate each inner rectangle point
        end
    end

    # Create nodes and assign them as free
    nodes = Node[Node(position, :free) for position in positions]

    # Connect nodes to form elements
    elements = connect_nodes(nodes, connectivity)

    # Assign IDs to nodes and elements
    for node in nodes
        node.id = :slab
    end
    for element in elements
        element.id = :beam
    end

    # Create the model and preprocess it
    model = Asap.Model(nodes, elements, Asap.AbstractLoad[])
    model = preprocess_model(model, i_columns=i_columns)

    # Plot the model if requested
    if plot
        plot_model(model)
    end

    return model

end
