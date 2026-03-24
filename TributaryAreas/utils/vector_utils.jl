"""
    crossproduct(a::Vector{Float64}, b::Vector{Float64}, o::Vector{Float64})

Calculate the cross product between vectors `ao` and `ob`, where `o` is the origin point.
"""
function crossproduct(a::Vector{Float64}, b::Vector{Float64}, o::Vector{Float64})
    x, y, x1, x2, y1, y2 = o[1], o[2], a[1], b[1], a[2], b[2]
    return (x2-x1) * (y-y1) - (x-x1) * (y2-y1)
end

"""
    crossproduct(a::Node, b::Node, o::Node)

Calculate the cross product between vectors `ao` and `ob`, where `o` is the origin point.
"""
function crossproduct(a::Node, b::Node, o::Node)
    x, y, x1, x2, y1, y2 = o.position[1], o.position[2], a.position[1], b.position[1], a.position[2], b.position[2]
    return crossproduct([x1, y1], [x2, y2], [x, y])
end

"""
    crossproduct(x1::Float64, x2::Float64, y1::Float64, y2::Float64)

Calculate the cross product of two 2D vectors given by their components.
"""
function crossproduct(x1::Float64, x2::Float64, y1::Float64, y2::Float64)
    return (x1 * y2) - (y1 * x2)
end

"""
    crossproduct(a::Vector{Float64}, b::Vector{Float64})

Calculate the cross product between two 2D vectors `a` and `b`.
"""
function crossproduct(a::Vector{Float64}, b::Vector{Float64})
    return crossproduct(a[1], b[1], a[2], b[2])
end

"""
    dotproduct(a::Vector{Float64}, b::Vector{Float64})

Calculate the dot product of two 2D vectors `a` and `b`.
"""
function dotproduct(a::Vector{Float64}, b::Vector{Float64})
    return a[1] * b[1] + a[2] * b[2]
end

"""
    project_vector(u::Vector{Float64}, v::Vector{Float64})

Project vector `u` onto vector `v` and return the projection and the perpendicular component.
"""
function project_vector(u::Vector{Float64}, v::Vector{Float64})
    projection = (dotproduct(u, v) / (v[1]^2 + v[2]^2)) * v
    perpendicular = u - projection
    return projection, perpendicular
end

"""
    unit_vector(x1::Float64, x2::Float64, y1::Float64, y2::Float64)

Calculate the unit vector from point `(x1, y1)` to point `(x2, y2)`.
"""
function unit_vector(x1::Float64, x2::Float64, y1::Float64, y2::Float64)
    vector = [(x2-x1), (y2-y1)]
    magnitude = sqrt((x2-x1)^2 + (y2-y1)^2)
    return vector / magnitude
end

"""
    unit_vector(vector::Vector{Float64})

Calculate the unit vector of a given 2D vector.
"""
function unit_vector(vector::Vector{Float64})
    magnitude = sqrt(vector[1]^2 + vector[2]^2)
    return vector / magnitude
end
