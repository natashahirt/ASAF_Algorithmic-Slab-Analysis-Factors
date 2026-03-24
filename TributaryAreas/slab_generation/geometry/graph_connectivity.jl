"""
    get_half_edges(nodes::Vector{Node}, elements::Vector{<:Element}; clockwise=true, return_hull=false)

Compute half-edges for a given set of nodes and elements. Optionally, return the hull.

# Arguments
- `nodes::Vector{Node}`: A vector of nodes.
- `elements::Vector{<:Element}`: A vector of elements.
- `clockwise::Bool`: If true, ensures edges are oriented clockwise.
- `return_hull::Bool`: If true, returns the hull.

# Returns
- `adjacency_dict::Dict`: A dictionary of node adjacencies.
- `half_edges::Set`: A set of half-edges.
- `nodes::Vector{Node}`: The input nodes.
- `i_hull_reversed::Vector` (optional): The reversed hull indices.
"""
function get_half_edges(nodes::Vector{Node}, elements::Vector{<:Element}; clockwise=true, return_hull=false)
    connectivity_matrix = build_connectivity_matrix(nodes, elements)
    i_hull, columns_hull = andrew_hull(nodes)
    hull_forward_edges, hull_backward_edges = determine_hull_edges(i_hull, clockwise)
    update_connectivity_matrix!(connectivity_matrix, nodes, elements, hull_forward_edges)
    adjacency_dict = build_adjacency_dict(connectivity_matrix, nodes, clockwise)
    half_edges = compute_half_edges(adjacency_dict, hull_forward_edges, hull_backward_edges)

    return return_hull ? (adjacency_dict, half_edges, nodes, reverse(i_hull)) : (adjacency_dict, half_edges, nodes)
end

"""
    build_connectivity_matrix(nodes, elements)

Constructs the initial connectivity matrix based on the given elements.
"""
function build_connectivity_matrix(nodes, elements)
    connectivity_matrix = zeros(Int64, lastindex(nodes), lastindex(nodes))
    for element in elements[:beam]
        connectivity_matrix[element.nodeStart.nodeID, element.nodeEnd.nodeID] = -1
        connectivity_matrix[element.nodeEnd.nodeID, element.nodeStart.nodeID] = 1
    end
    return connectivity_matrix
end

"""
    determine_hull_edges(i_hull, clockwise)

Determines the forward and backward hull edges based on the orientation.
"""
function determine_hull_edges(i_hull, clockwise)
    i_hull_reversed = reverse(i_hull)
    if clockwise
        hull_backward_edges = Set([[i_hull[i], i_hull[mod1(i+1, lastindex(i_hull))]] for i in 1:lastindex(i_hull)])
        hull_forward_edges = Set([[i_hull_reversed[i], i_hull_reversed[mod1(i+1, lastindex(i_hull))]] for i in 1:lastindex(i_hull)])
    else
        hull_forward_edges = Set([[i_hull[i], i_hull[mod1(i+1, lastindex(i_hull))]] for i in 1:lastindex(i_hull)])
        hull_backward_edges = Set([[i_hull_reversed[i], i_hull_reversed[mod1(i+1, lastindex(i_hull))]] for i in 1:lastindex(i_hull)])
    end
    return hull_forward_edges, hull_backward_edges
end

"""
    update_connectivity_matrix!(connectivity_matrix, nodes, elements, hull_forward_edges)

Updates the connectivity matrix with the hull edges.
"""
function update_connectivity_matrix!(connectivity_matrix, nodes, elements, hull_forward_edges)
    for edge in hull_forward_edges
        """if connectivity_matrix[edge[1], edge[2]] == 0 || connectivity_matrix[edge[2], edge[1]] == 0
            new_elements = connect_nodes(nodes, [edge])
            for element in new_elements; element.id = :beam; end
            append!(elements, new_elements)
        end"""
        connectivity_matrix[edge[1], edge[2]] = -1
        connectivity_matrix[edge[2], edge[1]] = 1
    end
end

"""
    build_adjacency_dict(connectivity_matrix, nodes, clockwise)

Builds and sorts the adjacency dictionary from the connectivity matrix.
"""
function build_adjacency_dict(connectivity_matrix, nodes, clockwise)
    adjacency_dict = Dict()
    for irow in axes(connectivity_matrix, 1)
        icol = findall(connectivity_matrix[irow, :] .!= 0)
        adjacency_dict[irow] = icol
    end

    for (i_node, i_adjacencies) in adjacency_dict
        if length(i_adjacencies) < 2 continue end
        node = nodes[i_node]
        adjacent_nodes = [nodes[i_adjacent] for i_adjacent in i_adjacencies]
        angles = [angle_pointpoint(node, adjacent_node, clockwise=clockwise) for adjacent_node in adjacent_nodes]
        perm = sortperm(angles)
        adjacency_dict[i_node] = i_adjacencies[perm]
    end
    return adjacency_dict
end

"""
    compute_half_edges(adjacency_dict, hull_forward_edges, hull_backward_edges)

Computes the set of half-edges from the adjacency dictionary and hull edges.
"""
function compute_half_edges(adjacency_dict, hull_forward_edges, hull_backward_edges)
    interior_half_edges = Set()
    for (i_node, i_adjacencies) in adjacency_dict
        for i_adjacency in i_adjacencies
            push!(interior_half_edges, [i_node, i_adjacency])
        end
    end
    interior_half_edges = setdiff(interior_half_edges, hull_backward_edges)
    return union(interior_half_edges, hull_forward_edges)
end