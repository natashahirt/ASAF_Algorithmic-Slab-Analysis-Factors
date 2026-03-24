"""
    J_<section>(x1::Real , ... , xn::Real ) 

Returns the polar moment of inertia of a section given n geometric properties for circles.
Torsional constant for non-circles.
Currently, different functions for each type of section

http://dir.cisc-icca.ca/files/technical/techdocs/updates/torsionprop.pdf
"""

function J_rect(h::Real , w::Real ) 
    # rectangular section
    # h = height, w = width

    if h == w # square

        J = 9/64 * h^4

    else # rectangle

        a = maximum([h,w])
        b = minimum([h,w])

        J = a * b ^3 / 3 - 0.21 * b^4 + 0.0175 * b^8 / a^4

    end

    return J

end

function J_I_symm(h::Real , w::Real , tw::Real , tf::Real ) 
    # I/H section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    J = J_I_asymm(h, w, w, tw, tf, tf)

    return J

end

function J_circ(D::Real) 
    # circular section
    # D = outer diameter

    J = π/32 * D^4

    return J

end

function J_circ_hollow(D::Real , d::Real ) 
    # hollow circular tube section
    # D = outer diameter, d = inner diameter

    J = π/32 * (D^4 - d^4)

    return J

end

function J_rect_hollow(H::Real , h::Real , W::Real , w::Real ) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, H = outer height, h = inner height

    t_h = (H-h)/2 # t
    t_w = (W-w)/2 # t1

    J = (2 * t_h * t_w * (H - t_h)^2 * (W - t_w)^2)/(H * t_h + W * t_w - t_h^2 - t_w^2)

    return J
    
end

function J_rect_hollow(H::Real , W::Real , t::Real ) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, t = wall thickness

    h = H - 2*t
    w = W - 2*t
    J = J_rect_hollow(H, h, W, w)

    return J
    
end

function J_U(h::Real , w::Real , tw::Real , tf::Real ) 
    # U section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    J = J_I_asymm(h, w, w, tw, tf, tf)

    return J
    
end

function J_T(h::Real , w::Real , tw::Real , tf::Real ) 
    # T section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    J = J_I_asymm(h, w, 0, tw, tf, 0)

    return J

end

function J_I_asymm(h::Real , w_t::Real , w_b::Real , tw::Real , tf_t::Real , tf_b::Real ) 
    # asymmetrical I/H section
    # h = height, w_t = width top, w_b = width bottom, tw = web thickness, tf_t = flange thickness top, tf_b = flange thickness bottom

    J = (((h - tf_t - tf_b) * tw^3) + (w_t * tf_t^3) + (w_b * tf_b^3))/3

    return J

end
