"""
    Iy_<section>(x1::Real , ... , xn::Real ) 

Returns the moment of inertia about a section's weak axis given n geometric properties
Currently, different functions for each type of section
"""

function Iy_rect(h::Real , w::Real ) 
    # rectangular section
    # h = height, w = width

    Iy = parallel_axis_I(h, w, 0.)

    return Iy

end

function Iy_I_symm(h::Real , w::Real , tw::Real , tf::Real ) 
    # I/H section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Iy = Iy_I_asymm(h, w, w, tw, tf, tf)

    return Iy

end

function Iy_circ(D::Real ) 
    # circular section
    # D = outer diameter

    Iy = D^4 * pi / 64 # strong axis

    return Iy

end

function Iy_circ_hollow(D::Real , d::Real ) 
    # hollow circular tube section
    # D = outer diameter, d = inner diameter

    Iy = Iy_circ(D) - Iy_circ(d)

    return Iy

end

function Iy_rect_hollow(H::Real , h::Real , W::Real , w::Real ) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, H = outer height, h = inner height

    Iy = Iy_rect(H,W) - Iy_rect(h,w)

    return Iy
    
end

function Iy_rect_hollow(H::Real , W::Real , t::Real ) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, t = wall thickness

    Iy = Iy_rect_hollow(H, (H - 2 * t), W, (W - 2 * t))

    return Iy

end

function Iy_U(h::Real , w::Real , tw::Real , tf::Real ) 
    # U section
    # h = height, w = width, tw = web thickness, tf = flange thickness

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

    # moments of inertia
    I_f = parallel_axis_I(tf, w, (x_f - x_c))
    I_w = parallel_axis_I(h_w, tw, (x_w - x_c))

    Iy = 2 * I_f + I_w

    return Iy
    
end

function Iy_T(h::Real , w::Real , tw::Real , tf::Real ) 
    # T section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Iy = Iy_I_asymm(h, w, 0., tw, tf, 0.)

    return Iy

end

@inline function Iy_I_asymm(h::Real , w_t::Real , w_b::Real , tw::Real , tf_t::Real , tf_b::Real ) 
    # asymmetrical I/H section
    # h = height, w_t = width top, w_b = width bottom, tw = web thickness, tf_t = flange thickness top, tf_b = flange thickness bottom

    h_w = h - tf_t - tf_b

    # moments of inertia
    I_b = parallel_axis_I(tf_b, w_b, 0.)
    I_w = parallel_axis_I(h_w, tw, 0.)
    I_t = parallel_axis_I(tf_t, w_t, 0.)

    Iy = I_b + I_w + I_t

    return Iy

end