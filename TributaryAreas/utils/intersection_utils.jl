"""
    line_line_intersect(a1::Float64, c1::Float64, a2::Float64, c2::Float64) -> Vector{Float64}

Calculate the intersection point of two lines given in the form `ax - by + c = 0`.

# Arguments
- `a1`, `c1`: Coefficients for the first line.
- `a2`, `c2`: Coefficients for the second line.

# Returns
- A vector `[x_intersect, y_intersect]` representing the intersection point.
"""
function line_line_intersect(a1::Float64, c1::Float64, a2::Float64, c2::Float64)
    x_intersect = (-c2 + c1) / (-a1 + a2)
    y_intersect = (c1 * a2 - a1 * c2) / (-a1 + a2)
    return [x_intersect, y_intersect]
end

"""
    ray_line_intersect(p::Vector{Float64}, r::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64}; param=0.5) -> Tuple{Vector{Float64}, Float64}

Calculate the intersection of a ray and a line segment.

# Arguments
- `p`: Starting point of the ray.
- `r`: Direction vector of the ray.
- `a`, `b`: Start and end points of the line segment.
- `param`: Parameter to determine the point along the ray.

# Returns
- A tuple `[point, distance]` where `point` is the intersection point and `distance` is the distance from `p` to `point`.
"""
function ray_line_intersect(p::Vector{Float64}, r::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64}; param=0.5)
    if a[2] == b[2] # Vertical line
        y = a[2]
        t = (a[2] - p[2]) / r[2]
        x = p[1] + t * r[1]
        intersect = [x, y]
    elseif a[1] == b[1] # Horizontal line
        x = a[1]
        t = (a[1] - p[1]) / r[1]
        y = p[2] + t * r[2]
        intersect = [x, y]
    else
        s = b - a
        q = a
        intersect, t = ray_ray_intersect(p, r, q, s)
    end

    if distance_pointpoint(a, b) < distance_pointpoint(intersect, a) || distance_pointpoint(a, b) < distance_pointpoint(intersect, b)
        return [nothing, 0.0]
    end

    point = p + t * param * r
    distance = distance_pointpoint(p, point)
    return [point, distance]
end

"""
    ray_line_intersect(p::Vector{Float64}, r::Vector{Float64}, edge::Vector{Vector{Float64}}; param=0.5) -> Tuple{Vector{Float64}, Float64}

Calculate the intersection of a ray and a line segment defined by an edge.

# Arguments
- `p`: Starting point of the ray.
- `r`: Direction vector of the ray.
- `edge`: A vector containing two vectors representing the start and end points of the line segment.
- `param`: Parameter to determine the point along the ray.

# Returns
- A tuple `[point, distance]` where `point` is the intersection point and `distance` is the distance from `p` to `point`.
"""
function ray_line_intersect(p::Vector{Float64}, r::Vector{Float64}, edge::Vector{Vector{Float64}}; param=0.5)
    a, b = edge
    return ray_line_intersect(p, r, a, b, param=param)
end

"""
    ray_ray_intersect(p::Vector{Float64}, r::Vector{Float64}, q::Vector{Float64}, s::Vector{Float64}) -> Tuple{Vector{Float64}, Float64}

Calculate the intersection of two rays.

# Arguments
- `p`: Starting point of the first ray.
- `r`: Direction vector of the first ray.
- `q`: Starting point of the second ray.
- `s`: Direction vector of the second ray.

# Returns
- A tuple `[intersection, t]` where `intersection` is the intersection point and `t` is the parameter for the first ray.
"""
function ray_ray_intersect(p::Vector{Float64}, r::Vector{Float64}, q::Vector{Float64}, s::Vector{Float64})
    t = crossproduct((q - p), s) / crossproduct(r, s)
    u = crossproduct((q - p), r) / crossproduct(r, s)
    x, y = p + t * r
    return [x, y], t
end

"""
    projects_to_line(line_start::Vector{Float64}, line_end::Vector{Float64}, point::Vector{Float64}) -> Bool

Check if a point projects onto a line segment.
"""
function project_to_line(line_start::Vector{Float64}, line_end::Vector{Float64}, point::Vector{Float64},)
    line_vector = line_end - line_start
    point_vector = point - line_start
    return project_to_line(line_vector, point_vector)
end

function project_to_line(line_vector::Vector{Float64}, point_vector::Vector{Float64})
    t = dot(point_vector, line_vector) / dot(line_vector, line_vector)
    return t
end

function project_to_line(line_start::Node, line_end::Node, point::Vector{Float64})
    return project_to_line(line_start.position[1:2], line_end.position[1:2], point)
end