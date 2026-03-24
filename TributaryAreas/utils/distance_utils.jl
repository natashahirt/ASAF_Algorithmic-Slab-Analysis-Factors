"""
    distance_pointpoint(x1::Float64, x2::Float64, y1::Float64, y2::Float64) -> Float64

Calculate the Euclidean distance between two points (x1, y1) and (x2, y2).
"""
function distance_pointpoint(x1::Float64, x2::Float64, y1::Float64, y2::Float64)
    return sqrt((x2 - x1)^2 + (y2 - y1)^2)
end

"""
    distance_pointpoint(a::Vector{Float64}, b::Vector{Float64}) -> Float64

Calculate the Euclidean distance between two points represented as vectors `a` and `b`.
"""
function distance_pointpoint(a::Vector{Float64}, b::Vector{Float64})
    return distance_pointpoint(a[1], b[1], a[2], b[2])
end

"""
    distance_pointpoint(a::Node, b::Node) -> Float64

Calculate the Euclidean distance between two `Node` objects `a` and `b`.
"""
function distance_pointpoint(a::Node, b::Node)
    x1, x2, y1, y2 = a.position[1], b.position[1], a.position[2], b.position[2]
    return distance_pointpoint(x1, x2, y1, y2)
end

"""
    get_length(edge::Vector{Vector{Float64}}) -> Float64

Calculate the length of an edge represented by a vector of two points.
"""
function get_length(edge::Vector{Vector{Float64}})
    x1, y1, x2, y2 = edge[1][1], edge[1][2], edge[2][1], edge[2][2]
    return distance_pointpoint(x1, x2, y1, y2)
end

"""
    get_length(element::T) where {T<:Element} -> Float64

Calculate the length of an element defined by its start and end nodes.
"""
function get_length(element::T) where {T<:Element}
    return distance_pointpoint(element.nodeStart, element.nodeEnd)
end

"""
    distance_pointline(x::Float64, y::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64) -> Float64

Calculate the perpendicular distance from a point (x, y) to a line defined by two points (x1, y1) and (x2, y2).
"""
function distance_pointline(x::Float64, y::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64)
    # Use the line equation to calculate the perpendicular distance
    distance = abs((x2 - x1) * (y1 - y) - (x1 - x) * (y2 - y1)) / sqrt((x2 - x1)^2 + (y2 - y1)^2)
    return distance
end

"""
    distance_pointline(a::Node, b::Node, point::Node) -> Float64

Calculate the perpendicular distance from a `Node` point to a line defined by two `Node` objects `a` and `b`.
"""
function distance_pointline(a::Node, b::Node, point::Node)
    x, y, x1, x2, y1, y2 = point.position[1], point.position[2], a.position[1], b.position[1], a.position[2], b.position[2]
    distance = distance_pointline(x, y, x1, x2, y1, y2)
    return distance
end

function distance_pointline(a::Node, b::Node, point::Vector{Float64})
    x, y, x1, x2, y1, y2 = point[1], point[2], a.position[1], b.position[1], a.position[2], b.position[2]
    distance = distance_pointline(x, y, x1, x2, y1, y2)
    return distance
end

"""
    get_distance_list(startpoints::Vector{Node}) -> Vector{Float64}

Calculate the distances between consecutive points in a list of `Node` objects, returning a list of distances.
"""
function get_distance_list(startpoints::Vector{Node})
    distance_list = Float64[]
    for i in 1:lastindex(startpoints)
        # Calculate distance between current point and the next, wrapping around at the end
        distance = distance_pointpoint(startpoints[i], startpoints[mod1(i + 1, lastindex(startpoints))])
        push!(distance_list, distance)
    end
    return distance_list
end
