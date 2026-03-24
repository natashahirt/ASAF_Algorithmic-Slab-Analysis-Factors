"""
    get_mean_node(nodes::Vector{Node}) -> Node

Calculate the mean position of a list of nodes and return a new Node
at that mean position with a fixed z-coordinate.
"""
function get_mean_node(nodes::Vector{Node})
    vectors = [[node.position[1], node.position[2]] for node in nodes]
    x_mean, y_mean = get_mean_point(vectors)
    mean_node = Node([x_mean, y_mean, 0], :zfixed)
    return mean_node
end

"""
    get_mean_point(nodes::Vector{Vector{T}}) where T <: Real -> Vector{Float64}

Calculate the mean x and y coordinates from a list of 2D points.
"""
function get_mean_point(nodes::Union{Vector{Vector{T}}, Vector{Vector}}) where T <: Real
    x_mean = round(mean([node[1] for node in nodes]), digits=2)
    y_mean = round(mean([node[2] for node in nodes]), digits=2)
    return [x_mean, y_mean]
end

"""
    closestpoint_pointline(x::Float64, y::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64) -> Vector{Float64}

Find the closest point on a line segment defined by (x1, y1) and (x2, y2) to a point (x, y).
"""
function closestpoint_pointline(x::Float64, y::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64)
    ab = [x2 - x1, y2 - y1]
    ap = [x1 - x, y1 - y]
    dot_product = ab[1] * ap[1] + ab[2] * ap[2]
    magnitude_ab = ab[1]^2 + ab[2]^2
    t = -dot_product / magnitude_ab

    if t < 0
        return [x1, y1]
    elseif t > 1
        return [x2, y2]
    end

    A = [x2 - x1 y2 - y1; y1 - y2 x2 - x1]
    B = -[(-x * (x2 - x1) - y * (y2 - y1)); (-y1 * (x2 - x1) + x1 * (y2 - y1))]
    x_closest, y_closest = A \ B
    return [x_closest, y_closest]
end

"""
    closestpoint_pointline(point::Node, nodeStart::Node, nodeEnd::Node) -> Vector{Float64}

Find the closest point on a line segment defined by nodeStart and nodeEnd to a point.
"""
function closestpoint_pointline(point::Node, nodeStart::Node, nodeEnd::Node)
    x, y = point.position[1], point.position[2]
    x1, y1 = nodeStart.position[1], nodeStart.position[2]
    x2, y2 = nodeEnd.position[1], nodeEnd.position[2]
    return closestpoint_pointline(x, y, x1, x2, y1, y2)
end

"""
    closestpoint_pointline(point::Vector{Float64}, nodeStart::Node, nodeEnd::Node) -> Vector{Float64}

Find the closest point on a line segment defined by nodeStart and nodeEnd to a 2D point.
"""
function closestpoint_pointline(point::Vector{Float64}, nodeStart::Node, nodeEnd::Node)
    x, y = point[1], point[2]
    x1, y1 = nodeStart.position[1], nodeStart.position[2]
    x2, y2 = nodeEnd.position[1], nodeEnd.position[2]
    return closestpoint_pointline(x, y, x1, x2, y1, y2)
end

"""
    param_pointline(x::Float64, y::Float64, nodeStart::Node, nodeEnd::Node) -> Float64

Calculate the parameter t for the projection of a point (x, y) onto a line segment
defined by nodeStart and nodeEnd.
"""
function param_pointline(x::Float64, y::Float64, nodeStart::Node, nodeEnd::Node)
    x1, y1 = nodeStart.position[1], nodeStart.position[2]
    x2, y2 = nodeEnd.position[1], nodeEnd.position[2]
    line_length = distance_pointpoint(x1, x2, y1, y2)
    distance_to_start = distance_pointpoint(x, x1, y, y1)
    param = distance_to_start / line_length
    return param
end

"""
    point_paramline(param::Float64, line::Vector{Vector{Float64}}) -> Vector{Float64}

Calculate the point on a line segment at a given parameter t.
"""
function point_paramline(param::Float64, line::Vector{Vector{Float64}})
    x1, x2, y1, y2 = line[1][1], line[2][1], line[1][2], line[2][2]
    x = (x2 - x1) * param + x1
    y = (y2 - y1) * param + y1
    return [x, y]
end

