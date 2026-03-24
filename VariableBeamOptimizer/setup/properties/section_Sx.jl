"""
    Sx_<section>(x1::Real, ... , xn::Real) 

Returns the elastic section modulus about a section's strong axis given n geometric properties
Sx = Ix / Y, where Y is the distance of the extreme fibre to the centroid

Currently, different functions for each type of section
"""

function Sx_rect(h::Real, w::Real) 
    # rectangular section
    # h = height, w = width

    Ix = Ix_rect(h,w)
    Y = h/2

    Sx = Ix/Y

    return Sx

end

function Sx_I_symm(h::Real, w::Real, tw::Real, tf::Real) 
    # I/H section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Ix = Ix_I_symm(h, w, tw, tf)
    Y = h/2

    Sx = Ix/Y

    return Sx

end

function Sx_circ(D::Real) 
    # circular section
    # D = outer diameter

    Ix = Ix_circ(D)
    Y = D/2

    Sx = Ix/Y

    return Sx

end

function Sx_circ_hollow(D::Real, d::Real) 
    # hollow circular tube section
    # D = outer diameter, d = inner diameter

    Ix = Ix_circ_hollow(D,d)
    Y = D/2

    Sx = Ix/Y

    return Sx

end

function Sx_rect_hollow(H::Real, h::Real, W::Real, w::Real) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, H = outer height, h = inner height

    Ix = Ix_rect_hollow(H,h,W,w)
    Y = H/2

    Sx = Ix/Y

    return Sx
    
end

function Sx_rect_hollow(H::Real, W::Real, t::Real) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, t = wall thickness

    Sx = Sx_rect_hollow(H,(H-2*t),W,(W-2*t))

    return Sx

end

function Sx_U(h::Real, w::Real, tw::Real, tf::Real) 
    # U section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Ix = Ix_U(h,w,tw,tf)
    Y = h/2

    Sx = Ix/Y

    return Sx
    
end

function Sx_T(h::Real, w::Real, tw::Real, tf::Real) 
    # T section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Ix = Ix_T(h,w,tw,tf)
    Y_b = elastic_neutral_axis_y(h, w, 0., tw, tf, 0.)

    Y_t = h - Y_b

    Sx_b = Ix / Y_b
    Sx_t = Ix / Y_t

    return Sx_t, Sx_b # top fibre, bottom fibre

end

function Sx_I_asymm(h::Real, w_t::Real, w_b::Real, tw::Real, tf_t::Real, tf_b::Real; verbose=false) 
    # asymmetrical I/H section
    # h = height, w_t = width top, w_b = width bottom, tw = web thickness, tf_t = flange thickness top, tf_b = flange thickness bottom

    Ix = Ix_I_asymm(h, w_t, w_b, tw, tf_t, tf_b)
    Y_b = elastic_neutral_axis_y(h, w_t, w_b, tw, tf_t, tf_b)
    Y_t = h - Y_b

    Sx_b = Ix / Y_b
    Sx_t = Ix / Y_t

    return Sx_t, Sx_b # top fibre, bottom fibre

end