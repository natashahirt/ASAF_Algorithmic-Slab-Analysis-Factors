"""
    check_collinearity(element_a::Element, element_b::Element) -> Bool

Check if two elements are collinear by comparing their edges.

# Arguments
- `element_a::Element`: The first element with start and end nodes.
- `element_b::Element`: The second element with start and end nodes.

# Returns
- `Bool`: `true` if the elements are collinear, `false` otherwise.
"""
function check_collinearity(element_a::Element, element_b::Element)
    # Extract the positions of the start and end nodes for both elements
    edge_a = [[element_a.nodeStart.position[1], element_a.nodeStart.position[2]], 
              [element_a.nodeEnd.position[1], element_a.nodeEnd.position[2]]]
    edge_b = [[element_b.nodeStart.position[1], element_b.nodeStart.position[2]], 
              [element_b.nodeEnd.position[1], element_b.nodeEnd.position[2]]]

    return check_collinearity(edge_a, edge_b)
end

"""
    check_collinearity(edge_a::Vector{Vector{Float64}}, edge_b::Vector{Vector{Float64}}; tol=1e-5) -> Bool

Check if two edges are collinear using their vector representation.

# Arguments
- `edge_a::Vector{Vector{Float64}}`: The first edge defined by two points.
- `edge_b::Vector{Vector{Float64}}`: The second edge defined by two points.
- `tol::Float64`: Tolerance for collinearity check (default is 1e-5).

# Returns
- `Bool`: `true` if the edges are collinear, `false` otherwise.
"""
function check_collinearity(edge_a::Vector{Vector{Float64}}, edge_b::Vector{Vector{Float64}}; tol=1e-5)
    # Calculate direction vectors for both edges
    vector_a = [edge_a[2][1] - edge_a[1][1], edge_a[2][2] - edge_a[1][2]]
    vector_b = [edge_b[2][1] - edge_b[1][1], edge_b[2][2] - edge_b[1][2]]
    
    # Calculate the cross product of the direction vectors
    cross_product = crossproduct(vector_a, vector_b)

    # Check if the cross product is within the tolerance
    return abs(cross_product) <= tol
end

"""
    check_collinearity(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}; tol=0.01) -> Bool

Check if three points are collinear by comparing the gradients of the lines they form.

# Arguments
- `a::Vector{Float64}`: The first point.
- `b::Vector{Float64}`: The second point.
- `c::Vector{Float64}`: The third point.
- `tol::Float64`: Tolerance for collinearity check (default is 0.01).

# Returns
- `Bool`: `true` if the points are collinear, `false` otherwise.
"""
function check_collinearity(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}; tol=0.01)
    # Calculate the gradients of the lines formed by the points
    m1 = line_gradient(a, b)
    m2 = line_gradient(b, c)
    
    # Check if the difference in gradients is within the tolerance
    return abs(m1 - m2) < tol
end
