"""
    get_side(a::Node, b::Node, point::Node) -> Int

Determine the relative position of a point with respect to a line segment
defined by nodes `a` and `b`. Returns -1 if the point is to the left of the line,
1 if to the right, and 0 if on the line.
"""
function get_side(a::Node, b::Node, point::Node)
    cross_product = crossproduct(a, b, point)
    return cross_product > 0 ? -1 : cross_product < 0 ? 1 : 0
end

"""
    split_sides(a::Node, b::Node, nodes::Vector{Node}) -> Tuple{Vector{Node}, Vector{Node}, Vector{Node}}

Splits the nodes into three groups: those on the left, those on the right, and those on the line
defined by nodes `a` and `b`.
"""
function split_sides(a::Node, b::Node, nodes::Vector{Node})
    left, right, on_line = Node[], Node[], Node[]
    for node in nodes
        side = get_side(a, b, node)
        if side == 1
            push!(right, node)
        elseif side == -1
            push!(left, node)
        else
            push!(on_line, node)  # keep nodes that are exactly on the line
        end
    end
    return left, right, on_line
end

"""
    get_hull(a::Node, b::Node, nodes::Vector{Node})

Recursively finds the convex hull points between nodes `a` and `b` from the given nodes.
"""
function get_hull(a::Node, b::Node, nodes::Vector{Node})
    if isempty(nodes)
        return
    end

    # Calculate distances of nodes to the line (a, b)
    distances = [distance_pointline(a, b, node) for node in nodes]
    
    # Select nodes that are off the line (positive distance) and capture those on the line
    off_line_nodes = [nodes[i] for i in eachindex(nodes) if distances[i] > 0]
    on_line_nodes = [nodes[i] for i in eachindex(nodes) if distances[i] == 0]

    # If no off-line nodes are found, return early
    if isempty(off_line_nodes)
        # Keep nodes on the line only if they haven't been added to the hull yet
        for node in on_line_nodes
            if node ∉ hull
                push!(hull, node)
            end
        end
        return
    end

    # Find the furthest point from the line
    p = argmax(distances)
    c = off_line_nodes[p]

    # Add the furthest point to the hull
    if c ∉ hull
        push!(hull, c)
    end

    # Split nodes into two halves for further recursion
    left1, right1, _ = split_sides(a, c, off_line_nodes)
    left2, right2, _ = split_sides(c, b, off_line_nodes)

    # Recursively process each half
    get_hull(a, c, right1)
    get_hull(c, b, right2)
end

"""
    quickhull_algorithm(nodes::Vector{Node}; clockwise=true) -> Tuple{Vector{Int}, Vector{Int}}

Computes the convex hull of a set of nodes using the QuickHull algorithm.
"""
function quickhull_algorithm(nodes::Vector{Node}; clockwise=true)
    @assert length(nodes) >= 3 "You don't have enough nodes to do a quickhull (less than 3)."

    # Sort nodes by x and y positions to find extreme points
    x_pos = [node.position[1] for node in nodes]
    y_pos = [node.position[2] for node in nodes]

    p_y = sortperm(y_pos)
    sorted_nodes = nodes[p_y]
    p_x = sortperm(x_pos)
    sorted_nodes = sorted_nodes[p_x]

    bl, tr = sorted_nodes[1], sorted_nodes[end]  # bottom-left and top-right points

    global hull = [bl, tr]

    # Trimmed nodes exclude the extreme points (bl and tr)
    trimmed_nodes = sorted_nodes[2:end-1]

    # Split into left and right sets of nodes relative to the line bl-tr
    left, right, _ = split_sides(bl, tr, trimmed_nodes)

    # Recursively find the hull points on both sides
    get_hull(bl, tr, right)
    get_hull(tr, bl, left)

    # Sort the final hull points by angle around the mean node to ensure correct order
    meannode = get_mean_node(hull)
    angles = angle_pointpoint(meannode, hull, clockwise=clockwise)
    p_sorted = sortperm(angles)
    hull .= hull[p_sorted]

    # Map hull points to their indices in the original node array
    i_hull = [hull_node.nodeID for hull_node in hull]
    columns_hull = copy(i_hull)

    # Optionally check for overlaps (if required by your problem)
    hull_edges = [[i_hull[i], i_hull[mod1(i+1, lastindex(i_hull))]] for i in 1:lastindex(i_hull)]
    for edge in hull_edges
        overlap_i = check_overlap(nodes, edge)
        append!(i_hull, overlap_i)
    end

    # Ensure no duplicate indices
    i_hull = collect(Set(i_hull))
    hull = [nodes[i] for i in i_hull]

    # Re-sort hull points after overlap handling
    meannode = get_mean_node(hull)
    angles = angle_pointpoint(meannode, hull, clockwise=clockwise)
    p_sorted = sortperm(angles)
    i_hull .= i_hull[p_sorted]

    return i_hull, columns_hull
end

"""
    andrew_hull(nodes::Vector{Node}) -> Vector{Int}

Computes the convex hull of a set of nodes using Andrew's monotone chain algorithm.
"""
function andrew_hull(points::Vector{Vector{Float64}}; clockwise=true)
    nodes = [Node([point[1], point[2], 0], :free) for (i, point) in enumerate(points)]
    for (i, node) in enumerate(nodes)
        node.nodeID = i
    end
    return andrew_hull(nodes, clockwise=clockwise)
end

function andrew_hull(x::Vector{Float64}, y::Vector{Float64}; clockwise=true)
    nodes = Node[Node([x[i], y[i], 0], :free) for i in eachindex(x)]
    for (i, node) in enumerate(nodes)
        node.nodeID = i
    end
    return andrew_hull(nodes, clockwise=clockwise)
end

function andrew_hull(nodes::Vector{Node}; clockwise=true)
    @assert length(nodes) >= 3 "You don't have enough nodes to compute a convex hull (less than 3)."

    # Sort nodes by x-coordinate (and y-coordinate as a tie breaker)
    sorted_nodes = sort(nodes, by = x -> (x.position[1], x.position[2]))

    # Function to check if three points make a clockwise, counterclockwise, or collinear turn
    function cross(o::Node, a::Node, b::Node)
        (a.position[1] - o.position[1]) * (b.position[2] - o.position[2]) -
        (a.position[2] - o.position[2]) * (b.position[1] - o.position[1])
    end

    # Build the lower hull
    lower_hull = Node[]
    for node in sorted_nodes
        while length(lower_hull) >= 2 && cross(lower_hull[end-1], lower_hull[end], node) <= 0
            pop!(lower_hull)  # Remove the last point if it's not a "left turn"
        end
        push!(lower_hull, node)
    end

    # Build the upper hull
    upper_hull = Node[]
    for node in reverse(sorted_nodes)
        while length(upper_hull) >= 2 && cross(upper_hull[end-1], upper_hull[end], node) <= 0
            pop!(upper_hull)  # Remove the last point if it's not a "left turn"
        end
        push!(upper_hull, node)
    end

    # Remove the last point of each half because it's repeated
    pop!(lower_hull)
    pop!(upper_hull)

    # Concatenate the lower and upper hulls to get the full convex hull
    full_hull = vcat(lower_hull, upper_hull)

    # Return the indices of the hull points in the original node array
    i_hull = [hull_node.nodeID for hull_node in full_hull]
    columns_hull = copy(i_hull)

    # Optionally check for overlaps (if required by your problem)
    hull_edges = [[i_hull[i], i_hull[mod1(i+1, lastindex(i_hull))]] for i in 1:lastindex(i_hull)]
    for edge in hull_edges
        overlap_i = check_overlap(nodes, edge)
        append!(i_hull, overlap_i)
    end

    # Ensure no duplicate indices
    i_hull = collect(Set(i_hull))
    hull = [nodes[i] for i in i_hull]

    # Re-sort hull points after overlap handling
    meannode = get_mean_node(hull)
    angles = angle_pointpoint(meannode, hull, clockwise=clockwise)
    p_sorted = sortperm(angles)
    i_hull .= i_hull[p_sorted]

    return i_hull, columns_hull
end