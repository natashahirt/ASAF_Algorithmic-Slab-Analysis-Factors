"""
    get_next_cycle_node(adjacency_dict, half_edges, i_visited)

Extracts the next node in a cycle from the adjacency dictionary.

# Arguments
- `adjacency_dict::Dict`: A dictionary representing node adjacencies.
- `half_edges::Set`: A set of half edges.
- `i_visited::Vector`: A vector of visited nodes.

# Returns
- The next node in the cycle or `nothing` if no valid cycle is found.
"""
function get_next_cycle_node(adjacency_dict::Dict, half_edges::Set, i_visited::Vector)
    i_current = i_visited[end]
    i_adjacencies = adjacency_dict[i_current]

    if length(i_visited) > 1
        # Adjust adjacency order to prioritize unvisited nodes
        i_circshift = findfirst(==(i_visited[end-1]), i_adjacencies)
        i_adjacencies = circshift(i_adjacencies, -i_circshift)
    end

    for i_adjacent in i_adjacencies
        if i_adjacent == i_visited[1]
            if is_valid_half_edge([i_current, i_adjacent], half_edges)
                return i_adjacent
            end
        elseif !(i_adjacent in i_visited)
            if is_valid_half_edge([i_current, i_adjacent], half_edges)
                return i_adjacent
            end
        end
    end

    return nothing
end

"""
    is_valid_half_edge(half_edge, half_edges)

Checks if a half edge is valid and removes it from the set if valid.

# Arguments
- `half_edge::Vector`: A vector representing the half edge.
- `half_edges::Set`: A set of half edges.

# Returns
- `true` if the half edge is valid, otherwise `false`.
"""
function is_valid_half_edge(half_edge::Vector, half_edges::Set)
    if half_edge in half_edges
        delete!(half_edges, half_edge)
        return true
    end
    return false
end

"""
    get_cycles(adjacency_dict, half_edges)

Finds all counter-clockwise cycles in the adjacency dictionary.

# Arguments
- `adjacency_dict::Dict`: A dictionary representing node adjacencies.
- `half_edges::Set`: A set of half edges.

# Returns
- A vector of paths representing the cycles found.
"""
function get_cycles(adjacency_dict::Dict, half_edges::Set)
    paths_found = Vector[]

    for start_key in keys(adjacency_dict)
        i_visited = [start_key]

        while length(i_visited) <= length(keys(adjacency_dict))
            i_next = get_next_cycle_node(adjacency_dict, half_edges, i_visited)

            if isnothing(i_next)
                break
            elseif i_next == i_visited[1]
                push!(paths_found, i_visited)
                i_visited = Int64[]
            end

            push!(i_visited, i_next)
        end
    end

    return paths_found
end

"""
    filter_valid_cycles(self, cycles)

Filter out invalid cycles based on node collinearity.

# Arguments
- `self`: The parameters for slab analysis.
- `cycles`: A vector of cycles, where each cycle is a vector of node indices.

# Returns
- `valid_cycles`: A vector of valid cycles where not all nodes are collinear.

This function iterates over each cycle and checks the collinearity of each triplet of consecutive nodes. A cycle is considered valid if at least one triplet of nodes is not collinear.
"""
function filter_valid_cycles(self, cycles)
    valid_cycles = Vector[]
    model = self.model

    for cycle in cycles
        # Skip cycles with less than 3 nodes
        length(cycle) < 3 && continue
        
        # Check if all edges have associated elements
        cycle_elements = get_cycle_elements(cycle, model.elements[:beam])
        isnothing(cycle_elements) && continue

        # Check that all nodes in the cycle are not hole nodes
        is_hole_cycle = check_for_hole_cycle(cycle, self.i_holes)

        # Skip this cycle it is a hole
        is_hole_cycle && continue

        slab_nodes = [model.nodes[i].position[1:2] for i in cycle]
        is_valid = false

        # Check each triplet of nodes
        for i in 1:lastindex(slab_nodes)
            node_a = slab_nodes[mod1(i-1, length(slab_nodes))]
            node_b = slab_nodes[i]
            node_c = slab_nodes[mod1(i+1, length(slab_nodes))]

            # Calculate area of triangle formed by three points
            # If area is non-zero, points are not collinear
            area = abs(
                (node_b[1] - node_a[1]) * (node_c[2] - node_a[2]) -
                (node_c[1] - node_a[1]) * (node_b[2] - node_a[2])
            )

            if area > 1e-10  # Use small tolerance to account for floating point errors
                is_valid = true
                break
            end
        end

        if is_valid
            push!(valid_cycles, cycle)
        end
    end

    return valid_cycles
end

function check_for_hole_cycle(cycle::Vector{Int64}, i_holes::Vector{Vector{Int64}})
    # Check if all nodes in any hole vector are present in the cycle
    for hole_indices in i_holes
        # If all indices in this hole vector are present in the cycle, it's a hole
        if all(i -> i in cycle, hole_indices)
            return true
        end
    end

    return false

end

