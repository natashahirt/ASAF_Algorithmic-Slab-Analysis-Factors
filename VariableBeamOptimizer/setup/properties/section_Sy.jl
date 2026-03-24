"""
    Sy_<section>(x1::Real, ... , xn::Real) 

Returns the elastic section modulus about a section's weak axis given n geometric properties
Sy = Ix / Y, where Y is the distance of the extreme fibre to the centroid

Currently, different functions for each type of section
"""

function Sy_rect(h::Real, w::Real) 
    # rectangular section
    # h = height, w = width

    Iy = Iy_rect(h,w)
    X = w / 2

    Sy = Iy / X

    return Sy

end

function Sy_I_symm(h::Real, w::Real, tw::Real, tf::Real) 
    # I/H section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Iy = Iy_I_symm(h, w, tw, tf)
    X = w / 2
 
    Sy = Iy / X

    return Sy

end

function Sy_circ(D::Real) 
    # circular section
    # D = outer diameter

    Iy = Iy_circ(D)
    X = D / 2

    Sy = Iy / X

    return Sy

end

function Sy_circ_hollow(D::Real, d::Real) 
    # hollow circular tube section
    # D = outer diameter, d = inner diameter

    Iy = Iy_circ_hollow
    X = D / 2

    Sy = Iy / X

    return Sy

end

function Sy_rect_hollow(H::Real, h::Real, W::Real, w::Real) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, H = outer height, h = inner height

    Iy = Iy_rect_hollow(H,h,W,w)
    X = W / 2

    Sy = Iy / X

    return Sy
    
end

function Sy_rect_hollow(H::Real, W::Real, t::Real) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, t = wall thickness

    Sy = Sy_rect_hollow(H,(H-2*t),W,(W-2*t))

    return Sy

end

function Sy_U(h::Real, w::Real, tw::Real, tf::Real) 
    # U section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Iy = Iy_U(h,w,tw,tf)
    
    X_l = elastic_neutral_axis_x(h, w, tw, tf)
    X_r = w - X_l

    Sy_l = Iy / X_l
    Sy_r = Iy / X_r

    return Sy_l, Sy_r
    
end

function Sy_T(h::Real, w::Real, tw::Real, tf::Real) 
    # T section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Iy = Iy_T(h,w,tw,tf)
    X = w / 2

    Sy = Iy / X

    return Sy

end

function Sy_I_asymm(h::Real, w_t::Real, w_b::Real, tw::Real, tf_t::Real, tf_b::Real; verbose=false) 
    # asymmetrical I/H section
    # h = height, w_t = width top, w_b = width bottom, tw = web thickness, tf_t = flange thickness top, tf_b = flange thickness bottom

    Iy = Iy_I_asymm(h, w_t, w_b, tw, tf_t, tf_b)
    X = max(w_t, w_b) / 2

    Sy = Iy / X

    return Sy

end