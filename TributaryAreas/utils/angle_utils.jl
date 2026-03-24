"""
    angle_pointpoint(x1::Float64, x2::Float64, y1::Float64, y2::Float64; clockwise=true, degrees=false)

Calculate the angle between two points in a 2D plane.

# Arguments
- `x1::Float64`: x-coordinate of the first point
- `x2::Float64`: x-coordinate of the second point
- `y1::Float64`: y-coordinate of the first point
- `y2::Float64`: y-coordinate of the second point
- `clockwise::Bool=true`: if true, angle is measured clockwise; otherwise, counterclockwise
- `degrees::Bool=false`: if true, return angle in degrees; otherwise, in radians

# Returns
- `angle::Float64`: The calculated angle
"""
function angle_pointpoint(x1::Float64, x2::Float64, y1::Float64, y2::Float64; clockwise=true, degrees=false)
    # Calculate the differences in x and y coordinates
    x_prime, y_prime = x2 - x1, y2 - y1
    
    # Calculate the angle (counterclockwise by default)
    angle = atan(y_prime, x_prime)
 
    # Adjust the angle if clockwise measurement is requested
    if clockwise
        angle = angle < 0 ? abs(angle) : 2pi - angle
    end

    # Convert to degrees if requested
    if degrees
        angle *= 180 / pi
    end
    
    return angle
end

"""
    angle_pointpoint(a::Vector{Float64}, b::Vector{Float64}; clockwise=true, degrees=false)

Calculate the angle between two points represented as vectors in a 2D plane.

# Arguments
- `a::Vector{Float64}`: First point as a vector [x, y]
- `b::Vector{Float64}`: Second point as a vector [x, y]
- `clockwise::Bool=true`: if true, angle is measured clockwise; otherwise, counterclockwise
- `degrees::Bool=false`: if true, return angle in degrees; otherwise, in radians

# Returns
- `angle::Float64`: The calculated angle
"""
function angle_pointpoint(a::Vector{Float64}, b::Vector{Float64}; clockwise=true, degrees=false)
    # Directly pass vector elements to the main function
    return angle_pointpoint(a[1], b[1], a[2], b[2], clockwise=clockwise, degrees=degrees)
end

"""
    angle_pointpoint(a::Node, b::Node; clockwise=true, degrees=false)

Calculate the angle between two nodes in a 2D plane.

# Arguments
- `a::Node`: First node with position attribute
- `b::Node`: Second node with position attribute
- `clockwise::Bool=true`: if true, angle is measured clockwise; otherwise, counterclockwise
- `degrees::Bool=false`: if true, return angle in degrees; otherwise, in radians

# Returns
- `angle::Float64`: The calculated angle
"""
function angle_pointpoint(a::Node, b::Node; clockwise=true, degrees=false)
    # Directly pass node position elements to the main function
    return angle_pointpoint(a.position[1], b.position[1], a.position[2], b.position[2], clockwise=clockwise, degrees=degrees)
end

"""
    angle_pointpoint(a::Node, nodes::Vector{Node}; clockwise=true, degrees=false)

Calculate the angles between a node and a list of nodes in a 2D plane.

# Arguments
- `a::Node`: Reference node with position attribute
- `nodes::Vector{Node}`: List of nodes to calculate angles with
- `clockwise::Bool=true`: if true, angle is measured clockwise; otherwise, counterclockwise
- `degrees::Bool=false`: if true, return angle in degrees; otherwise, in radians

# Returns
- `angles::Vector{Float64}`: The calculated angles
"""
function angle_pointpoint(a::Node, nodes::Vector{Node}; clockwise=true, degrees=false)
    # Preallocate the angles vector for efficiency
    angles = Vector{Float64}(undef, length(nodes))
    # Calculate angle for each node in the list
    for (i, node) in enumerate(nodes)
        angles[i] = angle_pointpoint(a, node, clockwise=clockwise, degrees=degrees)
    end
    return angles
end

"""
    angle_pointpoint(a::Vector{Float64}, points::Vector{Vector}; clockwise=true, degrees=false)

Calculate the angles between a point and a list of points in a 2D plane.

# Arguments
- `a::Vector{Float64}`: Reference point as a vector [x, y]
- `points::Vector{Vector}`: List of points to calculate angles with
- `clockwise::Bool=true`: if true, angle is measured clockwise; otherwise, counterclockwise
- `degrees::Bool=false`: if true, return angle in degrees; otherwise, in radians

# Returns
- `angles::Vector{Float64}`: The calculated angles
"""
function angle_pointpoint(a::Vector{Float64}, points::Vector{Vector}; clockwise=true, degrees=false)
    # Preallocate the angles vector for efficiency
    angles = Vector{Float64}(undef, length(points))
    # Calculate angle for each point in the list
    for (i, point) in enumerate(points)
        angles[i] = angle_pointpoint(a, point, clockwise=clockwise, degrees=degrees)
    end
    return angles
end
