"""
    generate_regular_polygon(diameter::Real, n_sides::Int64; centroid=false, i_columns=Int64[])

Generates a regular polygon model with a specified number of sides and diameter.

# Arguments
- `diameter::Real`: The diameter of the polygon.
- `n_sides::Int64`: The number of sides of the polygon.
- `centroid::Bool`: Optional. If `true`, includes the centroid in the model. Default is `false`.
- `i_columns::Vector{Int64}`: Optional. A vector of column indices. Default is an empty vector.
- `plot::Bool`: Whether to plot the model. Defaults to false.

# Returns
- `model`: A preprocessed model of the polygon with nodes and elements.
"""
function generate_regular_polygon(diameter::Real, n_sides::Int64; centroid=false, i_columns=Int64[], plot=false)

    positions = Vector{Vector{Float64}}(undef, n_sides)

    angle = 2 * π / n_sides # radians

    for i in 0:(n_sides-1)
        point_angle = i * angle
        positions[i+1] = [diameter/2 * cos(point_angle), diameter/2 * sin(point_angle), 0.0]
    end

    connectivity = [[i, mod1(i+1, n_sides)] for i in 1:n_sides]

    if centroid
        push!(positions, [0.0, 0.0, 0.0])
        append!(connectivity, [n_sides+1, i] for i in 1:n_sides)
    end

    nodes = Node[Node(position, :free) for position in positions]
    elements = connect_nodes(nodes, connectivity)

    for node in nodes
        node.id = :slab
    end
    for element in elements
        element.id = :beam
    end

    model = Asap.Model(nodes, elements, Asap.AbstractLoad[])
    model = preprocess_model(model, i_columns=i_columns)

    if plot
        plot_model(model)
    end

    return model
end
