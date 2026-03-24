# geometric functions

"""
    parallel_axis_I(b::Real, t::Real, y::Real)

Get the moment of inertia, according to parallel axis theorem, for rectangles.

# Arguments
- b::Real : breadth of rectangle (perpendicular to y)
- t::Real : thickness of rectangle (parallel to y)
- y::Real : distance from neutral axis
"""
function parallel_axis_I(b, t, y)

    A = b * t
    I = b * t^3 / 12 + A * y^2

    return I

end

"""
    elastic_neutral_axis_y(h::Real, w_t::Real, w_b::Real, tw::Real, tf_t::Real, tf_b::Real)

Get the y-neutral axis of an asymmetric-I-beam-inspired-geometry.
The neutral axis is calculated with respect to the bottom edge.

# Arguments
- h::Real : height of beam
- w_t::Real : top flange width
- w_b::Real : bottom flange width
- tw::Real : web thickness
- tf_t::Real : top flange thickness
- tf_b::Real : bottom flange thickness
"""
function elastic_neutral_axis_y(h::Real, w_t::Real, w_b::Real, tw::Real, tf_t::Real, tf_b::Real) 

    h_w = h - tf_t - tf_b

    # areas
    A_b = w_b * tf_b # area bottom
    A_w = h_w * tw # area web
    A_t = w_t * tf_t # area top

    A = A_b + A_w + A_t

    # distances to bottom edge
    y_b = tf_b / 2 # from centroid of lower flange
    y_w = tf_b + h_w / 2 # from centroid of web
    y_t = (h - tf_t / 2) # from centroid of top flange

    # distance to neutral axis from bottom edge
    y_c = ((A_b * y_b) + (A_w * y_w) + (A_t) * y_t) / A

    return y_c

end

"""
    elastic_neutral_axis_x(h::Real, w::Real, tw::Real, tf::Real)

Get the x-neutral axis of an channel-inspired-geometry.
The neutral axis is calculated with respect to the left edge (also the edge where the web is).

# Arguments
- h::Real : height of beam
- w::Real : width of beam
- tw::Real : web thickness
- tf::Real : flange thickness
"""
function elastic_neutral_axis_x(h::Real, w::Real, tw::Real, tf::Real) 

    h_w = h - 2 * tf

    # areas
    A_f = w * tf # area flange
    A_w = h_w * tw # area web

    A = 2 * A_f + A_w

    # distances to left edge
    x_f = w / 2 # from centroid of flange
    x_w = tw / 2 # from centroid of web

    # distance to neutral axis from bottom edge
    x_c = (2 * (A_f * x_f) + (A_w * x_w)) / A

    return x_c

end

"""
    plastic_neutral_axis_y(h::Real, w_t::Real, w_b::Real, tw::Real, tf_t::Real, tf_b::Real)

Get the plastic neutral axis of an I-beam-inspired geometry.
The plastic neutral axis is the point of the beam where areas above and below are equal
(i.e. areas are unweighted by their distance from the centroid.) The assumption is that the
compression area (Ac) and tension area (At) under a uniform stress field will be equal.

Ac * fy = At * fy ⇔ Ac = At

# Arguments
- h::Real : height of beam
- w_t::Real : top flange width
- w_b::Real : bottom flange width
- tw::Real : web thickness
- tf_t::Real : top flange thickness
- tf_b::Real : bottom flange thickness
"""
function plastic_neutral_axis_y(h::Real, w_t::Real, w_b::Real, tw::Real, tf_t::Real, tf_b::Real) 

    h_w = h - tf_t - tf_b

    # areas
    A_b = w_b * tf_b # area bottom
    A_w = h_w * tw # area web
    A_t = w_t * tf_t # area top

    A = A_b + A_w + A_t

    # find y_pna from bottom flange
    if A_b >= A / 2 # pna is in the bottom flange

        # y distance from bottom edge taken up by A/2
        y_pna = (A / 2) / w_b

    elseif A_t >= A / 2 # pna is in the top flange

        # y distance from bottom edge taken up by A/2
        y_pna = h - ((A / 2) / w_t)

    else # pna is in the web

        y_pna = tf_b + (A - 2 * w_b * tf_b) / (2 * tw)

    end

    return y_pna

end

"""
    plastic_neutral_axis_x(h::Real, w::Real, tw::Real, tf::Real)

Get the x-plastic-neutral axis of an channel-inspired-geometry (weak axis).
The neutral axis is calculated with respect to the left edge (also the edge where the web is).

# Arguments
- h::Real : height of beam
- w::Real : width of beam
- tw::Real : web thickness
- tf::Real : flange thickness
"""
function plastic_neutral_axis_x(h::Real, w::Real, tw::Real, tf::Real) 
    # h = height, w = width, tw = web thickness, tf = flange thickness

    w_trimmed = w - tw

    # areas
    A_f = w_trimmed * tf # area flange (without the web)
    A_w = h * tw # area web (full height)

    A = 2 * A_f + A_w

    # get x_pna
    if A_w >= A / 2 # neutral axis is in the web

        x_pna = (A / 2) / h

    else

        x_pna = w - A / (4 * tf)

    end
    
    return x_pna

end

# shared functions

"""
    RG(x1::Real, ... , xn::Real) 

Returns the gradius of gyration of a section relative to an axis (works for x and y).
One function for all sections. Input processed values.

# Arguments
- A::Float64 : area of section
- I::Float64 : moment of inertia about the desired axis (Ix, Iy)
"""

function RG(A::Real, I::Real)

    if A <= 0
        A = 1e-6
    end

    if I <= 0
        I = 1e-6
    end

    @assert A > 0 "Area is zero!"
    @assert I > 0 "Moment of inertia is zero!"

    return sqrt(I / A)

end