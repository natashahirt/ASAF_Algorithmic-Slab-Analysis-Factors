# ==============================================================================
# AISC Design Guide 9 — Torsion for Doubly-Symmetric I-Shapes (W-shapes)
# ==============================================================================
#
# Stresses computed:
#   - Pure (St. Venant) torsional shear: τ_t = G·t·θ'     (Eq. 4.1)
#   - Warping shear stress:  τ_ws = E·Sw·|θ'''| / t        (Eq. 4.2a)
#   - Warping normal stress: σ_w  = E·Wno·θ''              (Eq. 4.3a)
#
# Design checks (LRFD, §4.7.1):
#   - Normal stress yielding: |σ_b| + |σ_w| ≤ φ·Fy        (Eq. 4.12)
#   - Shear stress yielding:  τ_b + τ_t + τ_ws ≤ φ·0.6·Fy (Eq. 4.13)
#   - Interaction: (fun/(φFy))² + (fuv/(φ·0.6Fy))² ≤ 1.0  (Eq. 4.16a)
#
# All functions use raw kip/inch units (no Unitful).
# ==============================================================================

# ==============================================================================
# DG9 Torsional Properties
# ==============================================================================

"""
    dg9_Wno(s::I_symm) -> Real

Normalized warping function at flange tips: Wno = bf · ho / 4.
"""
dg9_Wno(s::I_symm) = s.bf * s.ho / 4

"""
    dg9_Sw1(s::I_symm) -> Real

Warping statical moment at web-flange junction: Sw1 = tf · bf² · ho / 16.
"""
dg9_Sw1(s::I_symm) = s.tf * s.bf^2 * s.ho / 16

"""
    dg9_torsional_parameter(s::I_symm) -> Real

Torsional parameter 'a' = √(E·Cw / (G·J)) per DG9 Eq. 3.4.
"""
function dg9_torsional_parameter(s::I_symm)
    E, G = s.material.E, s.material.G
    return sqrt(E * s.Cw / (G * s.J))
end

# ==============================================================================
# Internal helper — extract torsion properties as a NamedTuple
# ==============================================================================

function _torsion_props(s::I_symm)
    Wno = dg9_Wno(s)
    Sw1 = dg9_Sw1(s)
    return (d=s.d, bf=s.bf, tw=s.tw, tf=s.tf, ho=s.ho,
            J=s.J, Cw=s.Cw, Ix=s.Ix, Sx=s.Sx,
            Wno=Wno, Sw1=Sw1)
end

# ==============================================================================
# Torsional Functions — Loading Cases (DG9 Appendix B)
# ==============================================================================

"""
    torsion_case3_derivatives(z, L, T, a, G, J)

θ and derivatives for concentrated torque T at midspan, pinned-pinned.
Returns (θ, θp, θpp, θppp).
"""
function torsion_case3_derivatives(z::Real, L::Real, T::Real,
                                    a::Real, G_val::Real, J_val::Real)
    α = L / (2 * a)
    ζ = z / a
    GJ = G_val * J_val
    half_TGJ = T / (2 * GJ)
    cosh_α = cosh(α)

    if z ≤ L / 2
        sinh_ζ = sinh(ζ)
        cosh_ζ = cosh(ζ)
        θ    = half_TGJ * (z - a * sinh_ζ / cosh_α)
        θp   = half_TGJ * (1.0 - cosh_ζ / cosh_α)
        θpp  = -half_TGJ / a * (sinh_ζ / cosh_α)
        θppp = -half_TGJ / a^2 * (cosh_ζ / cosh_α)
    else
        z2 = L - z
        ζ2 = z2 / a
        sinh_ζ2 = sinh(ζ2)
        cosh_ζ2 = cosh(ζ2)
        θ    = half_TGJ * (z2 - a * sinh_ζ2 / cosh_α)
        θp   = -(half_TGJ * (1.0 - cosh_ζ2 / cosh_α))
        θpp  = half_TGJ / a * (sinh_ζ2 / cosh_α)
        θppp = -half_TGJ / a^2 * (cosh_ζ2 / cosh_α)
    end

    return (θ=θ, θp=θp, θpp=θpp, θppp=θppp)
end

"""
    torsion_case1_derivatives(z, L, t_per_length, a, G, J)

θ and derivatives for uniform torque, pinned-pinned.
"""
function torsion_case1_derivatives(z::Real, L::Real, t_per_length::Real,
                                    a::Real, G_val::Real, J_val::Real)
    GJ = G_val * J_val
    tGJ = t_per_length / GJ
    α = L / (2 * a)
    ζ = z / a
    cosh_α = cosh(α)
    cosh_shifted = cosh(ζ - α)

    θ    = tGJ * (z * (L - z) / 2 + a^2 * (cosh_shifted / cosh_α - 1))
    θp   = tGJ * ((L - 2*z) / 2 + a * sinh(ζ - α) / cosh_α)
    θpp  = tGJ * (-1.0 + cosh_shifted / cosh_α)
    θppp = tGJ / a * sinh(ζ - α) / cosh_α

    return (θ=θ, θp=θp, θpp=θpp, θppp=θppp)
end

# ==============================================================================
# Torsional Stress Calculations
# ==============================================================================

"""
    torsional_stresses(E, G, tf, tw, d, Ix, Wno, Sw1, θp, θpp, θppp; Vu=0)

Compute torsional stresses. All inputs and outputs in consistent units (ksi).
"""
function torsional_stresses(E::Real, G::Real,
                            tf::Real, tw::Real, d::Real, Ix::Real,
                            Wno::Real, Sw1::Real,
                            θp::Real, θpp::Real, θppp::Real;
                            Vu::Real=0.0)
    τ_t_flange = G * tf * abs(θp)
    τ_t_web    = G * tw * abs(θp)
    σ_w = E * Wno * θpp
    τ_ws = E * Sw1 * abs(θppp) / tf
    τ_b_web = abs(Vu) / (d * tw)

    return (τ_t_flange=τ_t_flange, τ_t_web=τ_t_web,
            σ_w=σ_w, τ_ws=τ_ws, τ_b_web=τ_b_web)
end

# ==============================================================================
# Design Checks (DG9 §4.7.1 — LRFD)
# ==============================================================================

"""
    check_torsion_yielding(σ_b, σ_w, τ_b, τ_t, τ_ws, Fy; φ=0.90) -> NamedTuple

Yielding check under combined bending + torsion per DG9 §4.7.1 (LRFD).
"""
function check_torsion_yielding(σ_b::Real, σ_w::Real,
                                 τ_b::Real, τ_t::Real, τ_ws::Real,
                                 Fy::Real; φ::Real=0.90)
    f_un = abs(σ_b) + abs(σ_w)
    f_uv = abs(τ_b) + abs(τ_t) + abs(τ_ws)

    φFy  = φ * Fy
    φFvy = φ * 0.6 * Fy

    normal_ok = f_un ≤ φFy
    shear_ok  = f_uv ≤ φFvy

    ir = (f_un / φFy)^2 + (f_uv / φFvy)^2
    interaction_ok = ir ≤ 1.0

    return (f_un=f_un, f_uv=f_uv, φFy=φFy, φFvy=φFvy,
            normal_ok=normal_ok, shear_ok=shear_ok,
            interaction_ratio=ir, interaction_ok=interaction_ok,
            ok=normal_ok && shear_ok && interaction_ok)
end

# ==============================================================================
# Full Torsion Design Check
# ==============================================================================

"""
    design_w_torsion(s::I_symm, Tu, Vu, Mu, L;
                     load_type=:concentrated_midspan) -> NamedTuple

Complete torsion design check per AISC Design Guide 9.

All arguments are in raw kip/inch units (Tu in kip·in, Vu in kip,
Mu in kip·in, L in inches). Returns stresses in ksi.

# Keyword Arguments
- `load_type`: `:concentrated_midspan` (DG9 Case 3) or `:uniform` (DG9 Case 1)
"""
function design_w_torsion(s::I_symm, Tu::Real, Vu::Real, Mu::Real, L::Real;
                          load_type::Symbol=:concentrated_midspan)
    E  = s.material.E
    Fy = s.material.Fy
    G  = s.material.G

    p = _torsion_props(s)
    a = sqrt(E * p.Cw / (G * p.J))

    Tu_abs = abs(Tu)
    Vu_abs = abs(Vu)
    Mu_abs = abs(Mu)

    if load_type == :concentrated_midspan
        d_mid = torsion_case3_derivatives(L/2, L, Tu_abs, a, G, p.J)
        d_sup = torsion_case3_derivatives(0.0, L, Tu_abs, a, G, p.J)
    elseif load_type == :uniform
        t_per_in = Tu_abs / L
        d_mid = torsion_case1_derivatives(L/2, L, t_per_in, a, G, p.J)
        d_sup = torsion_case1_derivatives(0.0, L, t_per_in, a, G, p.J)
    else
        error("Unsupported load_type: $load_type")
    end

    σ_b = Mu_abs / p.Sx

    s_mid = torsional_stresses(E, G, p.tf, p.tw, p.d, p.Ix,
                               p.Wno, p.Sw1, d_mid.θp, d_mid.θpp, d_mid.θppp)
    s_sup = torsional_stresses(E, G, p.tf, p.tw, p.d, p.Ix,
                               p.Wno, p.Sw1, d_sup.θp, d_sup.θpp, d_sup.θppp;
                               Vu=Vu_abs)

    f_un_mid = abs(σ_b) + abs(s_mid.σ_w)
    f_uv_mid = s_mid.τ_t_flange + s_mid.τ_ws

    f_un_sup = abs(s_sup.σ_w)
    f_uv_sup = s_sup.τ_b_web + s_sup.τ_t_flange + s_sup.τ_ws

    chk_mid = check_torsion_yielding(σ_b, s_mid.σ_w, 0.0,
                                      s_mid.τ_t_flange, s_mid.τ_ws, Fy)
    chk_sup = check_torsion_yielding(0.0, s_sup.σ_w, s_sup.τ_b_web,
                                      s_sup.τ_t_flange, s_sup.τ_ws, Fy)

    return (
        a = a,
        Wno = p.Wno,
        Sw1 = p.Sw1,
        σ_b = σ_b,
        σ_w_midspan = s_mid.σ_w,
        τ_t_midspan = s_mid.τ_t_flange,
        τ_ws_midspan = s_mid.τ_ws,
        f_un_midspan = f_un_mid,
        f_uv_midspan = f_uv_mid,
        σ_w_support = s_sup.σ_w,
        τ_t_support = s_sup.τ_t_flange,
        τ_ws_support = s_sup.τ_ws,
        τ_b_support = s_sup.τ_b_web,
        f_un_support = f_un_sup,
        f_uv_support = f_uv_sup,
        check_midspan = chk_mid,
        check_support = chk_sup,
        ok = chk_mid.ok && chk_sup.ok,
        θ_max_rad = d_mid.θ,
    )
end
