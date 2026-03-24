"""
    generate_rozvany(h::Real; index::Int64=0, i_columns=[1,2,3,4], plot=true)

Generates a Rozvany slab geometry based on the given height `h` and `index`.
The function returns a model with nodes and elements, optionally plotting the geometry.

# Arguments
- `h::Real`: The height of the slab.
- `index::Int64=0`: Determines the specific geometry configuration.
- `i_columns`: Indices of columns to be included in the model.
- `plot::Bool`: Whether to plot the model. Defaults to false.

# Returns
- `model`: A preprocessed model with nodes and elements.
"""
function generate_rozvany(h::Real; index::Int64=0, i_columns=[1,2,3,4], plot=false)
    # Initial square positions with columns
    positions = [[0., 0., 0.], [0., h, 0.], [h, h, 0.], [h, 0., 0.]]

    # Define connectivity and additional positions based on index
    if index == 0 || index == 8
        connectivity = [[1,2],[2,3],[3,4],[4,1]]
    elseif index == 1
        new_positions = [[0., h/2, 0.], [h/2, h, 0.], [h, h/2, 0.], [h/2, 0., 0.]]
        append!(positions, new_positions)
        connectivity = [[1,5],[5,2],[2,6],[6,3],[3,7],[7,4],[4,8],[8,1],[5,6],[6,7],[7,8],[8,5]]
    elseif index == 2
        new_positions = [[h/4,h/2,0.], [h/2,3h/4,0.], [3h/4,h/2,0.], [h/2,h/4,0.]]
        append!(positions, new_positions)
        connectivity = [[1,2],[2,3],[3,4],[4,1],[1,8],[1,5],[2,5],[2,6],[3,6],[3,7],[4,7],[4,8],[5,6],[6,7],[7,8],[8,5]]
    elseif index == 3
        new_positions = [[h/2,h,0.], [h/2,0.,0.], [h/4, sqrt(2)*h/4, 0.], [h/4, h-sqrt(2)*h/4, 0.], [3h/4, sqrt(2)*h/4, 0.], [3h/4, h-sqrt(2)*h/4, 0.]]
        append!(positions, new_positions)
        connectivity = [[1,2],[2,5],[5,3],[3,4],[4,6],[6,1],[1,7],[2,8],[3,10],[4,9],[5,8],[5,10],[6,7],[6,9],[7,8],[8,10],[10,9],[9,7]]
    elseif index == 4
        new_positions = [[0.,h/(1+sqrt(2)),0.], [h/(1+sqrt(2)),h-h/(2+sqrt(2)),0.], [h-h/(2+sqrt(2)),h/(1+sqrt(2)),0.], [h/(1+sqrt(2)),0.,0.]]
        append!(positions, new_positions)
        connectivity = [[1,5],[5,2],[2,3],[3,4],[4,8],[8,1],[5,6],[6,7],[7,8],[8,5],[2,6],[3,6],[3,7],[4,7]]
    elseif index == 5
        new_positions = [[0., h*(sqrt(2)-1),0.], [0., h - h*(sqrt(2)-1),0.], [h*(sqrt(2)-1), h, 0.], [h/sqrt(2), h-h*(sqrt(2)-1), 0.], [h/sqrt(2), h*(sqrt(2)-1), 0.], [h*(sqrt(2)-1), 0., 0.]]
        append!(positions, new_positions)
        connectivity = [[1,5],[5,6],[6,2],[2,7],[7,3],[3,4],[4,10],[10,1],[6,7],[7,8],[8,9],[9,10],[10,5],[5,9],[6,8],[3,8],[4,9]]
    elseif index == 6
        new_positions = [[0., h/2, 0.], [h/2, h, 0.], [h, h/2, 0.]]
        append!(positions, new_positions)
        connectivity = [[1,5],[5,2],[2,6],[6,3],[3,7],[7,4],[4,1],[5,6],[6,7],[7,5]]
    elseif index == 7
        new_positions = [[0., h-h/sqrt(2), 0.], [h/2, h-h/(2*sqrt(2)), 0.], [h, h-h/sqrt(2), 0.]]
        append!(positions, new_positions)
        connectivity = [[1,5],[5,2],[2,3],[3,7],[7,4],[4,1],[5,6],[2,6],[3,6],[7,6],[7,5]]
    elseif index == 9
        new_positions = [[0., h-h/(1+sqrt(2)), 0.], [h-2*(h/(2+sqrt(2))), h, 0.], [h-h/(2+sqrt(2)), h-h/(1+sqrt(2)), 0.], [h-h/(2+sqrt(2)), 0., 0.]]
        append!(positions, new_positions)
        connectivity = [[1,5],[5,2],[2,6],[6,3],[3,4],[4,8],[8,1],[5,6],[6,7],[7,5],[3,7],[8,7]]
    elseif index == 10
        new_positions = [[h/4,h,0.], [3h/4,h,0.], [3h/4,0.,0.], [h/4,0.,0.]]
        append!(positions, new_positions)
        connectivity = [[1,2],[2,5],[5,6],[6,3],[3,4],[4,7],[7,8],[8,1],[5,8],[6,7]]
    elseif index == 11
        new_positions = [[h/4,h/2,0.], [h/2,3h/4,0.], [3h/4,h/2,0.], [3h/4,0.,0.], [h/4,0.,0.]]
        append!(positions, new_positions)
        connectivity = [[1,2],[2,3],[3,4],[4,8],[8,9],[9,1],[2,5],[2,6],[3,6],[3,7],[7,8],[9,5],[5,6],[6,7],[7,5]]
    elseif index == 12
        new_positions = [[h/4,h/2,0.], [h/2,3h/4,0.], [3h/4,h/2,0.], [3h/4,sqrt(2)*h/4,0.], [h/2,0.,0.], [h/4,sqrt(2)*h/4,0.]]
        append!(positions, new_positions)
        connectivity = [[1,2],[2,3],[3,4],[4,9],[9,1],[1,10],[2,5],[2,6],[3,6],[3,7],[4,8],[9,8],[9,10],[5,6],[6,7],[7,8],[8,10],[10,5],[5,7]]
    elseif index == 13
        new_positions = [[h/(2+sqrt(2)),0.,0.], [h/(2+sqrt(2)),h/(1+sqrt(2)),0.], [h-h/(1+sqrt(2)),h-h/(2+sqrt(2)),0.], [h,h/(1+sqrt(2)),0.]]
        append!(positions, new_positions)
        connectivity = [[1,2],[2,3],[3,8],[8,4],[4,5],[5,1],[5,6],[6,2],[2,7],[6,7],[3,7],[7,8],[8,6]]
    elseif index == 14
        new_positions = [[h/4, h-sqrt(2)*h/4, 0.], [h/2,h,0.], [3h/4, h-sqrt(2)*h/4,0.], [3h/4,0.,0.], [h/4,0.,0.]]
        append!(positions, new_positions)
        connectivity = [[1,2],[2,6],[6,3],[3,4],[4,8],[8,9],[9,1],[9,5],[2,5],[5,6],[6,7],[7,3],[7,8],[7,5]]
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

    if plot
        plot_model(model)
    end

    return model
end