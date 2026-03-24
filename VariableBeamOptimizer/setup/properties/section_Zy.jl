"""
    Zy_<section>(x1::Real, ... , xn::Real) 

Returns the plastic section modulus about a section's weak axis given n geometric properties
Zy = ∑A_i_left_of_pna * y_i_left_of_pna + ∑A_i_right_of_pna * y_i_right_of_pna, where pna is the plastic neutral axis

The idea being that, for a constant stress across the section
- A_compression * fy = A_tension * fy ⇔ A_compression = A_tension

This is only really relevant for steel sections, since other materials don't have the necessary
ductility to fail in such a uniform way.

Currently, different functions for each type of section
"""

function Zy_rect(h::Real, w::Real) 
    # rectangular section
    # h = height, w = width

    Zy = h * w^2 / 4

    return Zy

end

function Zy_I_symm(h::Real, w::Real, tw::Real, tf::Real) 
    # I/H section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Zy = Zy_I_asymm(h, w, w, tw, tf, tf)

    return Zy

end

function Zy_circ(D::Real) 
    # circular section
    # D = outer diameter

    Zy = 4 * (D/2)^3 / 3

    return Zy

end

function Zy_circ_hollow(D::Real, d::Real) 
    # hollow circular tube section
    # D = outer diameter, d = inner diameter

    Zy = Zy_circ(D) - Zy_circ(d)

    return Zy

end

function Zy_rect_hollow(H::Real, h::Real, W::Real, w::Real) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, H = outer height, h = inner height

    Zy = Zy_rect(H,W) - Zy_rect(h,w)

    return Zy
    
end

function Zy_rect_hollow(H::Real, W::Real, t::Real) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, t = wall thickness

    Zy = Zy_rect_hollow(H, (H - 2 * t), W, (W - 2 * t))

    return Zy

end

function Zy_U(h::Real, w::Real, tw::Real, tf::Real) 
    # U section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    w_trimmed = w - tw

    # areas
    A_f = w_trimmed * tf # area flange (without the web)
    A_w = h * tw # area web (full height)

    A = 2 * A_f + A_w

    # distances to left edge
    x_f = w_trimmed / 2 + tw # from centroid of flange
    x_w = tw / 2 # from centroid of web

    # get x_pna
    if A_w >= A / 2 # neutral axis is in the web

        x_pna = (A / 2) / h

        # region left of x_pna
        A_l_web = x_pna * h
        x_pna_l_web = x_pna / 2

        # region right of x_pna
        A_r_web = (tw - x_pna) * h
        x_pna_r_web = (tw - x_pna) / 2
        A_r_flange = w_trimmed * tf
        x_pna_r_flange = x_f - x_pna

        println("Ar web : $A_r_web, Ar flange: $A_r_flange")

        # get Zy
        Zy = A_l_web * x_pna_l_web + A_r_web * x_pna_r_web + 2 * (A_r_flange * x_pna_r_flange)

    else

        x_pna = w - A / (4 * tf)

        # region left of x_pna
        A_l_web = A_w
        x_pna_l_web = x_pna - x_w
        A_l_flange = (x_pna - tw) * tf
        x_pna_l_flange = (x_pna - tw) / 2

        # region right of x_pna
        A_r_flange = (w - x_pna) * tf
        x_pna_r_flange = (w - x_pna) / 2

        # get Zy
        Zy = A_l_web * x_pna_l_web + 2 * (A_l_flange * x_pna_l_flange) + 2 * (A_r_flange * x_pna_r_flange)

    end

    return Zy
    
end

function Zy_T(h::Real, w::Real, tw::Real, tf::Real) 
    # T section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Zy = Zy_I_asymm(h, w, 0, tw, tf, 0)

    return Zy

end

function Zy_I_asymm(h::Real, w_t::Real, w_b::Real, tw::Real, tf_t::Real, tf_b::Real) 
    # asymmetrical I/H section
    # h = height, w_t = width top, w_b = width bottom, tw = web thickness, tf_t = flange thickness top, tf_b = flange thickness bottom
    # note that I/H is only asymmetrical about the x axis, not the y axis
    
    h_w = h - tf_t - tf_b # height of web alone

    x_pna = h/2 # halfway 

    # region left of x_pna (region right of x_pna is identical)

    A_l_flange_t = w_t / 2 * tf_t
    x_pna_l_flange_t = w_t / 4
    A_l_web = h_w * tw / 2
    x_pna_l_web = tw / 4
    A_l_flange_b = w_b / 2 * tf_b
    x_pna_l_flange_b = w_b / 4

    # get Zy

    Zy = 2 * (A_l_flange_t * x_pna_l_flange_t + A_l_web * x_pna_l_web + A_l_flange_b * x_pna_l_flange_b)

    return Zy

end