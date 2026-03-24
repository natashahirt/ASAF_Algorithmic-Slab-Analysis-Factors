"""
    A_<section>(x1::Real , ... , xn::Real ) 

Returns the area of a section given n geometric properties
Currently, different functions for each type of section
"""

function A_rect(h::Real , w::Real ) 
    # rectangular section
    # h = height, w = width

    A = h * w

    return A

end

function A_I_symm(h::Real , w::Real , tw::Real , tf::Real ) 
    # I/H section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    A = A_I_asymm(h, w, w, tw, tf, tf)

    return A

end

function A_circ(D::Real ) 
    # circular section
    # D = outer diameter

    A = pi * (D/2)^2

    return A

end

function A_circ_hollow(D::Real , d::Real ) 
    # hollow circular tube section
    # D = outer diameter, d = inner diameter

    A = A_circ(D) - A_circ(d)

    return A

end

function A_rect_hollow(H::Real , h::Real , W::Real , w::Real ) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, H = outer height, h = inner height

    A = A_rect(H,W) - A_rect(h,w)

    return A
    
end

function A_rect_hollow(H::Real , W::Real , t::Real ) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, t = wall thickness

    A = A_rect_hollow(H, (H - 2 * t), W, (W - 2 * t))

    return A

end

function A_U(h::Real , w::Real , tw::Real , tf::Real ) 
    # U section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    A = A_I_asymm(h, w, w, tw, tf, tf)

    return A
    
end

function A_T(h::Real , w::Real , tw::Real , tf::Real ) 
    # T section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    A = A_I_asymm(h, w, 0., tw, tf, 0.)

    return A

end

function A_I_asymm(h::Real , w_t::Real , w_b::Real , tw::Real , tf_t::Real , tf_b::Real ) 
    # asymmetrical I/H section
    # h = height, w_t = width top, w_b = width bottom, tw = web thickness, tf_t = flange thickness top, tf_b = flange thickness bottom

    A_t = w_t * tf_t
    h_w = h - tf_t - tf_b
    A_w = h_w * tw 
    A_b = w_b * tf_b

    A = A_t + A_w + A_b

    return A

end
