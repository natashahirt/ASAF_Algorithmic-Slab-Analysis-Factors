"""
    line_gradient(a::Vector{Float64}, b::Vector{Float64})

Calculate the gradient (slope) of the line between two points `a` and `b`.
"""
function line_gradient(a::Vector{Float64}, b::Vector{Float64})
    delta_x = b[1] - a[1]
    delta_y = b[2] - a[2]
    return delta_y / delta_x
end

"""
    line_gradient(a::Node, b::Node)

Calculate the gradient (slope) of the line between two `Node` objects `a` and `b`.
"""
function line_gradient(a::Node, b::Node)
    a_pos = a.position
    b_pos = b.position
    return line_gradient([a_pos[1], a_pos[2]], [b_pos[1], b_pos[2]])
end

"""
    line_equation(a::Vector{Float64}, b::Vector{Float64})

Calculate the line equation parameters (slope `m` and y-intercept `c`) for the line between two points `a` and `b`.
Returns `m, c` where the line is in the form `y = mx + c`.
"""
function line_equation(a::Vector{Float64}, b::Vector{Float64})
    m = line_gradient(a, b)
    c = a[2] - a[1] * m  # y-intercept
    return m, c
end

"""
    line_equation(a::Node, b::Node)

Calculate the line equation parameters (slope `m` and y-intercept `c`) for the line between two `Node` objects `a` and `b`.
"""
function line_equation(a::Node, b::Node)
    a_pos = a.position
    b_pos = b.position
    return line_equation([a_pos[1], a_pos[2]], [b_pos[1], b_pos[2]])
end

"""
    perpendicularline_equation(element::Element, point::Vector{Float64})

Calculate the line equation parameters (slope `m_perp` and y-intercept `c_perp`) for the line perpendicular to the line defined by `element` and passing through `point`.
"""
function perpendicularline_equation(element::Element, point::Vector{Float64})
    a = element.nodeStart
    b = element.nodeEnd
    m_perp = -1 / line_gradient(a, b)
    c_perp = point[2] - m_perp * point[1]
    return m_perp, c_perp
end

"""
    perpendicularline_equation(a::Node, b::Node, point::Node)

Calculate the line equation parameters (slope `m_perp` and y-intercept `c_perp`) for the line perpendicular to the line between `Node` objects `a` and `b`, passing through `point`.
"""
function perpendicularline_equation(a::Node, b::Node, point::Node)
    m = line_gradient(a, b)
    m_perp = -1 / m
    c_perp = point.position[2] - m_perp * point.position[1]
    return m_perp, c_perp
end
