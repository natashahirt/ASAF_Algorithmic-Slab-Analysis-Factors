"""
    line_line_bisector(x::Float64, y::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64; param=0.5)

Calculate the bisector of two lines given their endpoints and a midpoint.

# Arguments
- `x`, `y`: Coordinates of the midpoint.
- `x1`, `y1`: Coordinates of the first endpoint.
- `x2`, `y2`: Coordinates of the second endpoint.
- `param`: Parameter to control the weight of the bisector direction.

# Returns
- A vector representing the bisector point.
"""
function line_line_bisector(x::Float64, y::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64; param=0.5)
    # Calculate unit vectors for each line
    unit_vector_1 = unit_vector(x, x1, y, y1) 
    unit_vector_2 = unit_vector(x, x2, y, y2)

    # Check if lines are collinear
    if check_collinearity([[x1,y1],[x,y]],[[x,y],[x2,y2]])
        # If collinear, return a perpendicular vector
        perp_vector = [unit_vector_1[2], -unit_vector_1[1]]
        bisector_point = [x, y] + perp_vector
        return bisector_point
    end

    # Calculate the bisector point
    bisector_point = unit_vector_1 * param + unit_vector_2 * (1-param) + [x, y]
    return bisector_point
end

"""
    line_line_bisector(o::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64}; param=0.5)

Overloaded function to calculate the bisector using vector inputs.

# Arguments
- `o`, `a`, `b`: Vectors representing the midpoint and endpoints.
- `param`: Parameter to control the weight of the bisector direction.

# Returns
- A vector representing the bisector point.
"""
function line_line_bisector(o::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64}; param=0.5)
    x, y, x1, y1, x2, y2 = o[1], o[2], a[1], a[2], b[1], b[2]
    return line_line_bisector(x, y, x1, x2, y1, y2, param=param)
end

"""
    line_line_bisector(o::Node, a::Node, b::Node; param=0.5)

Overloaded function to calculate the bisector using Node inputs.

# Arguments
- `o`, `a`, `b`: Nodes representing the midpoint and endpoints.
- `param`: Parameter to control the weight of the bisector direction.

# Returns
- A vector representing the bisector point.
"""
function line_line_bisector(o::Node, a::Node, b::Node; param=0.5)
    x, y, x1, y1, x2, y2 = o.position[1], o.position[2], a.position[1], a.position[2], b.position[1], b.position[2]
    return line_line_bisector(x, y, x1, x2, y1, y2, param=param)
end

"""
    get_bisector_intersections(bisectors::Vector{Vector}, startpoints::Vector{Vector})

Find intersections of bisectors with other bisectors.

# Arguments
- `bisectors`: A vector of bisector vectors.
- `startpoints`: A vector of starting points for each bisector.

# Returns
- A tuple containing the bisectors and their intersection points.
"""
function get_bisector_intersections(bisectors::Vector{Vector}, startpoints::Vector{Vector})
    intersections = Vector[]
    valid_bisectors = [bisector for bisector in bisectors if bisector != [0., 0.]]
    n_bisectors_to_check = Int(ceil(length(valid_bisectors) / 4))

    for i in 1:length(bisectors)
        intersection_points = Vector[]
        intersection_params = Float64[]

        if bisectors[i] == [0., 0.]
            push!(intersections, startpoints[i])
        else
            for j in (i-n_bisectors_to_check):(i+n_bisectors_to_check)
                j = mod1(j, length(bisectors))
                if j == i || bisectors[j] == [0., 0.]
                    continue
                end

                intersection_point, intersection_param = ray_ray_intersect(startpoints[i], bisectors[i], startpoints[j], bisectors[j])
                if intersection_param >= 0
                    push!(intersection_points, intersection_point)
                    push!(intersection_params, intersection_param)
                end
            end

            perm = sortperm(intersection_params)
            intersection_points = intersection_points[perm]

            if isempty(intersection_points)
                push!(intersections, startpoints[i])
                continue
            end

            push!(intersections, intersection_points[1])
        end
    end

    return bisectors, intersections
end

"""
    test_convexity(points::Vector{Vector})

Test if a polygon defined by a set of points is convex.

# Arguments
- `points`: A vector of points defining the polygon.

# Returns
- `true` if the polygon is convex, `false` otherwise.
"""
function test_convexity(points::Vector{Vector})
    crossproduct_z = []

    for (i, point) in enumerate(points)
        node_a = points[mod1(i-1, length(points))]
        node_o = point
        node_b = points[mod1(i+1, length(points))]
        push!(crossproduct_z, crossproduct(node_a, node_b, node_o))
    end

    # Use tolerance for floating point comparison
    tol = 1e-10
    pos_check = all(x -> x >= -tol, crossproduct_z)
    neg_check = all(x -> x <= tol, crossproduct_z)
    zero_check = all(x -> abs(x) <= tol, crossproduct_z)
    convex = pos_check || neg_check || zero_check

    return convex
end

"""
    test_convexity(nodes::Vector{Node})

Overloaded function to test convexity using Node inputs.

# Arguments
- `nodes`: A vector of Nodes defining the polygon.

# Returns
- `true` if the polygon is convex, `false` otherwise.
"""
function test_convexity(nodes::Vector{Node})
    points = [ [node.position[1], node.position[2]] for node in nodes ]
    return test_convexity(points)
end

"""
    rotate_point(centroid::Vector{Float64}, point::Vector{Float64}, theta_degrees::Float64)

Rotate a point around a centroid by a given angle.

# Arguments
- `centroid`: The centroid around which to rotate.
- `point`: The point to rotate.
- `theta_degrees`: The angle in degrees to rotate the point.

# Returns
- The rotated point.
"""
function rotate_point(centroid::Vector{Float64}, point::Vector{Float64}, theta_degrees::Float64)
    return_3d = false

    if length(point) == 3
        point = [point[1], point[2]]
        return_3d = true
    end

    radians = deg2rad(theta_degrees)
    centred_p = point - centroid
    rotation_matrix = [cos(radians) -sin(radians); sin(radians) cos(radians)]
    rotated_p = (rotation_matrix * centred_p) + centroid

    if return_3d
        rotated_p = [rotated_p[1], rotated_p[2], 0.]
    end

    return rotated_p
end

"""
    check_congruent_shape(distance_list::Vector, distance_lists::Vector)

Check if a shape is congruent to any in a list of shapes.

# Arguments
- `distance_list`: A vector of distances defining the shape.
- `distance_lists`: A vector of vectors, each defining a shape.

# Returns
- A tuple containing the index of the congruent shape, the indices of the matching points, and a boolean indicating if the shape is reversed.
"""
function check_congruent_shape(distance_list::Vector, distance_lists::Vector)
    reversed = false
    
    for i in 1:lastindex(distance_lists)
        if length(distance_list) != length(distance_lists[i])
            continue
        elseif sort(distance_list) != sort(distance_lists[i])
            continue
        else
            n = lastindex(distance_list)
            # Check for shifts
            for shift in 0:n-1
                shifted_list = circshift(distance_lists[i], shift)
                if shifted_list == distance_list
                    indices = circshift(collect(1:n), shift)
                    return i, indices, reversed
                end
            end

            # Check for reverse
            reversed_list = reverse(distance_lists[i])
            if reversed_list == distance_list
                reversed = true
                indices = reverse(collect(1:n))
                return i, indices, reversed
            end

            # Check for reverse + shift
            for shift in 0:n-1
                shifted_reversed_list = circshift(reversed_list, shift)
                if shifted_reversed_list == distance_list
                    reversed = true
                    indices = circshift(reverse(collect(1:n)), shift)
                    return i, indices, reversed
                end
            end
        end
    end

    return -1, [], reversed
end

"""
    get_polygon_area(nodes::Vector{Node})

Calculates the area of a polygon using the shoelace algorithm.

# Arguments
- `nodes::Vector{Node}`: A vector of `Node` objects representing the vertices of the polygon. Each `Node` should have a `position` attribute that provides the x and y coordinates.

# Returns
- `Float64`: The area of the polygon.

# Usage
This function is used to compute the area of a polygon defined by a sequence of nodes, typically in geometric computations.
"""
function get_polygon_area(nodes::Vector{Node})
    # Use the shoelace algorithm to find the polygon area

    n_vertices = lastindex(nodes)
    vertex_x = [node.position[1] for node in nodes]
    vertex_y = [node.position[2] for node in nodes]
    sum_1 = 0.
    sum_2 = 0.

    for i in 1:lastindex(nodes)
        i_plus = mod1(i+1, lastindex(nodes))
        sum_1 += vertex_x[i] * vertex_y[i_plus]
        sum_2 += vertex_y[i] * vertex_x[i_plus]
    end

    area = abs(sum_1 - sum_2) / 2

    return area
end
