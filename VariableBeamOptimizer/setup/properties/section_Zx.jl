"""
    Zx_<section>(x1::Real, ... , xn::Real) 

Returns the plastic section modulus about a section's strong axis given n geometric properties
Zx = ∑A_i_above_pna * y_i_above_pna + ∑A_i_below_pna * y_i_below_pna, where pna is the plastic neutral axis

The idea being that, for a constant stress across the section
- A_compression * fy = A_tension * fy ⇔ A_compression = A_tension

This is only really relevant for steel sections, since other materials don't have the necessary
ductility to fail in such a uniform way.

Currently, different functions for each type of section
"""

function Zx_rect(h::Real, w::Real) 
    # rectangular section
    # h = height, w = width

    Zx = w * h^2 / 4

    return Zx

end

function Zx_I_symm(h::Real, w::Real, tw::Real, tf::Real) 
    # I/H section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Zx = Zx_I_asymm(h, w, w, tw, tf, tf)

    return Zx

end

function Zx_circ(D::Real) 
    # circular section
    # D = outer diameter

    Zx = 4 * (D/2)^3 / 3

    return Zx

end

function Zx_circ_hollow(D::Real, d::Real) 
    # hollow circular tube section
    # D = outer diameter, d = inner diameter

    Zx = Zx_circ(D) - Zx_circ(d)

    return Zx

end

function Zx_rect_hollow(H::Real, h::Real, W::Real, w::Real) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, H = outer height, h = inner height

    Zx = Zx_rect(H,W) - Zx_rect(h,w)

    return Zx
    
end

function Zx_rect_hollow(H::Real, W::Real, t::Real) 
    # hollow rectangular tube section
    # W = outer width, w = inner width, t = wall thickness

    Zx = Zx_rect_hollow(H, (H - 2 * t), W, (W - 2 * t))

    return Zx

end

function Zx_U(h::Real, w::Real, tw::Real, tf::Real) 
    # U section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Zx = Zx_I_asymm(h, w, w, tw, tf, tf)

    return Zx
    
end

function Zx_T(h::Real, w::Real, tw::Real, tf::Real) 
    # T section
    # h = height, w = width, tw = web thickness, tf = flange thickness

    Zx = Zx_I_asymm(h, w, 0, tw, tf, 0)

    return Zx

end

function Zx_I_asymm(h::Real, w_t::Real, w_b::Real, tw::Real, tf_t::Real, tf_b::Real) 

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

    # find y_pna from bottom flange
    if A_b >= A / 2 # pna is in the bottom flange

        # y distance from bottom edge taken up by A/2
        y_pna = (A / 2) / w_b

        # region below y_pna
        A_b_bottom_flange = y_pna * w_b 
        y_pna_b_bottom_flange = y_pna / 2

        # region above y_pna (bottom flange, web)
        A_t_bottom_flange = (tf_b - y_pna) * w_b
        y_pna_t_bottom_flange = (tf_b - y_pna) / 2
        A_t_web = A_w
        y_pna_t_web = y_w - y_pna
        A_t_top_flange = A_t
        y_pna_t_top_flange = y_t - y_pna

        # get Zx
        Zx = A_b_bottom_flange * y_pna_b_bottom_flange + A_t_bottom_flange * y_pna_t_bottom_flange + A_t_web * y_pna_t_web + A_t_top_flange * y_pna_t_top_flange

    elseif A_t >= A / 2 # pna is in the top flange

        # y distance from bottom edge taken up by A/2
        y_pna = h - ((A / 2) / w_t)
        y_pna_from_top = (A / 2) / w_t

        # region below y_pna
        A_b_bottom_flange = A_b
        y_pna_b_bottom_flange = y_pna - y_b
        A_b_web = A_w
        y_pna_b_web = y_pna - y_w
        A_b_top_flange = abs(tf_t - y_pna_from_top) * w_t
        y_pna_b_top_flange = abs(tf_t - y_pna_from_top) / 2

        # region above y_pna
        A_t_top_flange = y_pna_from_top * w_t
        y_pna_t_top_flange = y_pna_from_top / 2

        # get Zx
        Zx = A_b_bottom_flange * y_pna_b_bottom_flange + A_b_web * y_pna_b_web + A_b_top_flange * y_pna_b_top_flange + A_t_top_flange * y_pna_t_top_flange

    else # pna is in the web

        y_pna = tf_b + (A - 2 * w_b * tf_b) / (2 * tw)

        # region below y_pna
        A_b_bottom_flange = A_b
        y_pna_b_bottom_flange = y_pna - y_b
        A_b_web = (y_pna - tf_b) * tw
        y_pna_b_web = (y_pna - tf_b) / 2

        # region above y_pna
        A_t_web = (h - tf_t - y_pna) * tw
        y_pna_t_web = (h - tf_t - y_pna) / 2
        A_t_top_flange = A_t
        y_pna_t_top_flange = y_t - y_pna

        # get Zx
        Zx = A_b_bottom_flange * y_pna_b_bottom_flange + A_b_web * y_pna_b_web + A_t_web * y_pna_t_web + A_t_top_flange * y_pna_t_top_flange

    end

    return Zx

end