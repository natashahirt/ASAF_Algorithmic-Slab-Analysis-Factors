"""
    generate_rectangle(h::Real, b::Real; i_columns=[1,2,3,4])

Generate a rectangular slab geometry with given height `h` and base `b`.
The optional parameter `i_columns` specifies the indices of columns.

# Arguments
- `h::Real`: Height of the rectangle.
- `b::Real`: Base of the rectangle.
- `i_columns`: Indices of columns (default is `[1,2,3,4]`).

# Returns
- `model`: A preprocessed model of the rectangle.
"""
function generate_rectangle(h::Real, b::Real; i_columns=[1,2,3,4], plot=false)
    # Define positions of the rectangle's corners
    positions = [[0., 0., 0.], [0., h, 0.], [b, h, 0.], [b, 0., 0.]]

    # Define connectivity between the corners
    connectivity = [[1,2], [2,3], [3,4], [4,1]]

    # Create nodes and elements
    nodes = Node[Node(position, :free) for position in positions]
    elements = connect_nodes(nodes, connectivity)

    # Assign IDs to nodes and elements
    for node in nodes node.id = :slab end
    for element in elements element.id = :beam end

    # Create and preprocess the model
    model = Asap.Model(nodes, elements, Asap.AbstractLoad[])
    model = preprocess_model(model, i_columns=i_columns)

    if plot
        plot_model(model)
    end

    return model
end

"""
    generate_rectangle(h::Real, b::Real, n::Int64; i_columns=[1,2,3,4])

Generate a rectangular slab geometry with girders, given height `h`, base `b`, 
and number of girders `n`. The optional parameter `i_columns` specifies the indices of columns.

# Arguments
- `h::Real`: Height of the rectangle.
- `b::Real`: Base of the rectangle.
- `n::Int64`: Number of girders.
- `i_columns`: Indices of columns (default is `[1,2,3,4]`).
- `plot::Bool`: Whether to plot the model. Defaults to false.

# Returns
- `model`: A preprocessed model of the rectangle with girders.
"""
function generate_rectangle(h::Real, b::Real, n::Int64; i_columns=[1,2,3,4], plot=false)
    # Define positions of the rectangle's corners
    corners = Vector[[0., 0., 0.], [0., h, 0.], [b, h, 0.], [b, 0., 0.]]

    # Calculate y-positions for girders
    girder_y = [i * h/(n+1) for i in 1:n]

    # Define positions for left and right girders
    left_girders = [[0., y, 0.] for y in girder_y]
    right_girders = [[b, y, 0.] for y in girder_y]

    # Indexing for girders
    left_girder_i = [i + 4 for i in 1:n]
    right_girder_i = [i + 4 + n for i in 1:n]

    # Combine all positions and indices
    i_s = [[1,2,3,4]; [left_girder_i; right_girder_i]]
    positions = Vector[corners; left_girders; right_girders]

    # Define connectivity for girders and edges
    girder_connectivities = [[left_girder_i[i], right_girder_i[i]] for i in 1:n]
    left_edge_connectivities = [[[1, left_girder_i[1]], [left_girder_i[end], 2]]; [[left_girder_i[i], left_girder_i[i+1]] for i in 1:n-1]]
    right_edge_connectivities = [[[4, right_girder_i[1]], [right_girder_i[end], 3]]; [[right_girder_i[i], right_girder_i[i+1]] for i in 1:n-1]]

    # Combine all connectivities
    connectivity = [[[2,3], [4,1]]; girder_connectivities; left_edge_connectivities; right_edge_connectivities]

    # Create nodes and elements
    nodes = Node[Node(position, :free) for position in positions]
    elements = connect_nodes(nodes, connectivity)

    # Assign IDs to nodes and elements
    for node in nodes node.id = :slab end
    for element in elements element.id = :beam end

    # Create and preprocess the model
    model = Asap.Model(nodes, elements, Asap.AbstractLoad[])
    model = preprocess_model(model, i_columns=i_columns)

    if plot
        plot_model(model)
    end

    return model
end
