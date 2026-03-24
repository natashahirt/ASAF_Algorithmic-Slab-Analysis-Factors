"""
    Ix_<section>(x1::Real , ... , xn::Real ) 

Returns the moment of inertia about a section's strong axis given n geometric properties
Currently, different functions for each type of section
"""

function Ix_rect(h::Real , w::Real ) 
    # rectangular section
    # h = height, w = width

    Ix = parallel_axis_I(w, h, 0)

    return Ix

end

function Ix_I_symm(h::Real , w::Real , tw::Real , tf::Real ) 
    # I/H section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Ix = Ix_I_asymm(h, w, w, tw, tf, tf)

    return Ix

end

function Ix_circ(D::Real ) 
    # circular section
    # D = outer diameter

    Ix = D^4 * pi / 64 # strong axis

    return Ix

end

function Ix_circ_hollow(D::Real , d::Real ) 
    # hollow circular tube section
    # D = outer diameter, d = inner diameter

    Ix = Ix_circ(D) - Ix_circ(d)

    return Ix

end

function Ix_rect_hollow(H::Real , h::Real , W::Real , w::Real ) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, H = outer height, h = inner height

    Ix = Ix_rect(H,W) - Ix_rect(h,w)

    return Ix
    
end

function Ix_rect_hollow(H::Real , W::Real , t::Real ) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, t = wall thickness

    Ix = Ix_rect_hollow(H, (H - 2 * t), W, (W - 2 * t))

    return Ix

end

function Ix_U(h::Real , w::Real , tw::Real , tf::Real ) 
    # U section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Ix = Ix_I_symm(h,w,tw,tf)

    return Ix
    
end

function Ix_T(h::Real , w::Real , tw::Real , tf::Real ) 
    # T section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Ix = Ix_I_asymm(h, w, 0., tw, tf, 0.)

    return Ix

end

function Ix_I_asymm(h::Real , w_t::Real , w_b::Real , tw::Real , tf_t::Real , tf_b::Real ) 
    # asymmetrical I/H section
    # h = height, w_t = width top, w_b = width bottom, tw = web thickness, tf_t = flange thickness top, tf_b = flange thickness bottom

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
    y_c = ((A_b * y_b) + (A_w * y_w) + (A_t * y_t)) / A

    # moments of inertia
    I_b = parallel_axis_I(w_b, tf_b, (y_b - y_c))
    I_w = parallel_axis_I(tw, h_w, (y_w - y_c))
    I_t = parallel_axis_I(w_t, tf_t, (y_t - y_c))

    Ix = I_b + I_w + I_t

    return Ix

end