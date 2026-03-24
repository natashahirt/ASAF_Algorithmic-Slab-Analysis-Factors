"""
    I_symm <: AbstractSection

Doubly-symmetric wide-flange I-section with AISC 360-16 design checks.

Input geometry uses AISC naming: `d` (depth), `bf` (flange width),
`tw` (web thickness), `tf` (flange thickness). All derived section
properties, slenderness classifications, and capacity values are
computed on construction.

# Fields (input)
- `d`   : Overall depth of section
- `bf`  : Flange width
- `tw`  : Web thickness
- `tf`  : Flange thickness
- `material` : `Metal` with E, G, Fy, Fu, ρ, ν, units
- `Lb`  : Unbraced length for LTB (0 = fully braced)

# Fields (computed) — see `get_properties` for full list
"""
mutable struct I_symm <: AbstractSection

    # --- Input geometry (AISC naming) ---
    d::Real         # overall depth
    bf::Real        # flange width
    tw::Real        # web thickness
    tf::Real        # flange thickness
    material::Metal
    Lb::Real        # unbraced length for LTB

    # --- Computed geometric properties ---
    A::Real         # cross-sectional area
    Ix::Real        # strong-axis moment of inertia
    Iy::Real        # weak-axis moment of inertia
    I_yc::Real      # compression-flange Iy about its own centroid
    J::Real         # St. Venant torsional constant
    Cw::Real        # warping constant
    Sx::Real        # strong-axis elastic section modulus
    Sy::Real        # weak-axis elastic section modulus
    Zx::Real        # strong-axis plastic section modulus
    Zy::Real        # weak-axis plastic section modulus
    rx::Real        # strong-axis radius of gyration
    ry::Real        # weak-axis radius of gyration
    ho::Real        # distance between flange centroids (d - tf)
    rts::Real       # effective radius of gyration for LTB (F2-7)

    # --- Neutral axes ---
    na_x::Real      # elastic NA, x-direction (= bf/2 for symmetric)
    na_y::Real      # elastic NA, y-direction (= d/2 for symmetric)
    pna_x::Real     # plastic NA, x-direction
    pna_y::Real     # plastic NA, y-direction

    # --- Slenderness (Table B4.1b, flexure) ---
    λ_f::Real       # flange width-thickness ratio bf/(2tf)
    λp_f::Real      # compact limit (flange)
    λr_f::Real      # noncompact limit (flange)
    λ_w::Real       # web width-thickness ratio h_w/tw
    λp_w::Real      # compact limit (web)
    λr_w::Real      # noncompact limit (web)
    slenderness_f::Symbol   # :compact, :noncompact, or :slender
    slenderness_w::Symbol

    # --- LTB parameters (Chapter F) ---
    Lp::Real        # limiting unbraced length for yielding (F2-5)
    Lr::Real        # limiting unbraced length for inelastic LTB (F2-6)
    F_cr::Real      # elastic LTB critical stress (F2-4)

    # --- Nominal capacities ---
    Mn::Real        # nominal strong-axis moment capacity
    Mn_weak::Real   # nominal weak-axis moment capacity (F6)
    Vn::Real        # nominal strong-axis shear capacity
    Vn_weak::Real   # nominal weak-axis shear capacity (G7)

    # --- Visualization ---
    flat_coords::Vector{Vector{Float64}}

    function I_symm(d::Real, bf::Real, tw::Real, tf::Real;
                    material::Metal=steel_ksi, Lb::Real=0)
        z = 0.0
        self = new(
            d, bf, tw, tf, material, Lb,
            z, z, z, z, z, z, z, z, z, z, z, z, z, z,  # A..rts (14)
            z, z, z, z,                                   # na_x..pna_y (4)
            z, z, z, z, z, z,                             # λ_f..λr_w (6)
            :compact, :compact,                           # slenderness_f, slenderness_w
            z, z, z,                                      # Lp, Lr, F_cr
            z, z, z, z,                                   # Mn, Mn_weak, Vn, Vn_weak
            Vector{Float64}[]                             # flat_coords
        )
        get_properties!(self)
        return self
    end
end

# Backward-compatible aliases: let downstream code use old field names
# `h` and `w` are the most commonly accessed old names.
Base.getproperty(s::I_symm, sym::Symbol) = begin
    if sym === :h
        getfield(s, :d)
    elseif sym === :w
        getfield(s, :bf)
    elseif sym === :RGx
        getfield(s, :rx)
    elseif sym === :RGy
        getfield(s, :ry)
    else
        getfield(s, sym)
    end
end

Base.setproperty!(s::I_symm, sym::Symbol, val) = begin
    if sym === :h
        setfield!(s, :d, val)
    elseif sym === :w
        setfield!(s, :bf, val)
    elseif sym === :RGx
        setfield!(s, :rx, val)
    elseif sym === :RGy
        setfield!(s, :ry, val)
    else
        setfield!(s, sym, val)
    end
end

# ============================================================================
# Property computation
# ============================================================================

"""
    get_properties!(self::I_symm)

Compute all derived section properties, slenderness classifications,
LTB parameters, and nominal capacities from the four input dimensions.
"""
function get_properties!(self::I_symm)
    d  = getfield(self, :d)
    bf = getfield(self, :bf)
    tw = getfield(self, :tw)
    tf = getfield(self, :tf)
    Lb = getfield(self, :Lb)
    mat = getfield(self, :material)

    # --- Cross-section properties ---
    A   = A_I_symm(d, bf, tw, tf)
    Ix  = Ix_I_symm(d, bf, tw, tf)
    Iy  = Iy_I_symm(d, bf, tw, tf)
    I_yc = parallel_axis_I(tf, bf, 0)
    J   = J_I_symm(d, bf, tw, tf)
    Sx  = Sx_I_symm(d, bf, tw, tf)
    Sy  = Sy_I_symm(d, bf, tw, tf)
    Zx  = Zx_I_symm(d, bf, tw, tf)
    Zy  = Zy_I_symm(d, bf, tw, tf)
    rx  = RG(A, Ix)
    ry  = RG(A, Iy)
    ho  = d - tf              # distance between flange centroids
    Cw  = Iy * ho^2 / 4      # warping constant

    rts_val = Sx > 0 ? sqrt(sqrt(Iy * Cw) / Sx) : 0.0

    na_x  = bf / 2
    na_y  = d / 2
    pna_x = bf / 2
    pna_y = d / 2

    setfield!(self, :A, A)
    setfield!(self, :Ix, Ix)
    setfield!(self, :Iy, Iy)
    setfield!(self, :I_yc, I_yc)
    setfield!(self, :J, J)
    setfield!(self, :Cw, Cw)
    setfield!(self, :Sx, Sx)
    setfield!(self, :Sy, Sy)
    setfield!(self, :Zx, Zx)
    setfield!(self, :Zy, Zy)
    setfield!(self, :rx, rx)
    setfield!(self, :ry, ry)
    setfield!(self, :ho, ho)
    setfield!(self, :rts, rts_val)
    setfield!(self, :na_x, na_x)
    setfield!(self, :na_y, na_y)
    setfield!(self, :pna_x, pna_x)
    setfield!(self, :pna_y, pna_y)

    # --- Slenderness classification (Table B4.1b) ---
    E, Fy = mat.E, mat.Fy
    h_w = d - 2 * tf

    λ_f  = bf / (2 * tf)
    λp_f = 0.38 * sqrt(E / Fy)
    λr_f = 1.0 * sqrt(E / Fy)
    class_f = λ_f > λr_f ? :slender : (λ_f > λp_f ? :noncompact : :compact)

    λ_w  = h_w / tw
    λp_w = 3.76 * sqrt(E / Fy)
    λr_w = 5.70 * sqrt(E / Fy)
    class_w = λ_w > λr_w ? :slender : (λ_w > λp_w ? :noncompact : :compact)

    setfield!(self, :λ_f, λ_f)
    setfield!(self, :λp_f, λp_f)
    setfield!(self, :λr_f, λr_f)
    setfield!(self, :λ_w, λ_w)
    setfield!(self, :λp_w, λp_w)
    setfield!(self, :λr_w, λr_w)
    setfield!(self, :slenderness_f, class_f)
    setfield!(self, :slenderness_w, class_w)

    # --- LTB parameters (Chapter F2) ---
    if Lb > 0
        Lp, Lr, F_cr = _compute_LTB_params(ry, J, Sx, ho, rts_val, Cw, Iy, Lb, E, Fy)
    else
        Lp = Lr = F_cr = 0.0
    end
    setfield!(self, :Lp, Lp)
    setfield!(self, :Lr, Lr)
    setfield!(self, :F_cr, F_cr)

    # --- Nominal capacities ---
    Mn      = _compute_Mn_strong(Zx, Sx, λ_f, λp_f, λr_f, λ_w, Lb, Lp, Lr, F_cr, E, Fy)
    Mn_weak = _compute_Mn_weak(Zy, Sy, λ_f, E, Fy)
    Vn      = _compute_Vn_strong(d, tw, tf, E, Fy)
    Vn_weak = _compute_Vn_weak(bf, tf, Fy)

    setfield!(self, :Mn, Mn)
    setfield!(self, :Mn_weak, Mn_weak)
    setfield!(self, :Vn, Vn)
    setfield!(self, :Vn_weak, Vn_weak)

    # --- Visualization ---
    setfield!(self, :flat_coords, _compute_coords(d, bf, tw, tf))
end

# Keep old name as alias for backward compatibility in callers
get_properties(self::I_symm) = get_properties!(self)

# ============================================================================
# LTB parameters (F2-5, F2-6, F2-4)
# ============================================================================

function _compute_LTB_params(ry, J, Sx, ho, rts, Cw, Iy, Lb, E, Fy)
    c = 1.0  # doubly symmetric
    Cb = 1.0 # conservative default

    Lp = 1.76 * ry * sqrt(E / Fy)

    jc_term = (J * c) / (Sx * ho)
    Lr = 1.95 * rts * (E / (0.7 * Fy)) *
         sqrt(jc_term + sqrt(jc_term^2 + 6.76 * (0.7 * Fy / E)^2))

    lb_rts = Lb / rts
    F_cr = Cb * π^2 * E / lb_rts^2 *
           sqrt(1 + 0.078 * jc_term * lb_rts^2)

    return Lp, Lr, F_cr
end

# ============================================================================
# Nominal flexural strength — Strong axis (F2/F3)
# ============================================================================

function _compute_Mn_strong(Zx, Sx, λ_f, λp_f, λr_f, λ_w, Lb, Lp, Lr, F_cr, E, Fy)
    Mp = Fy * Zx
    Mn = Mp

    # LTB (F2-2 / F2-3)
    if Lb > 0
        Cb = 1.0
        if Lb > Lr
            Mn = min(Mn, F_cr * Sx)
        elseif Lb > Lp
            Mn_LTB = Cb * (Mp - (Mp - 0.7 * Fy * Sx) * ((Lb - Lp) / (Lr - Lp)))
            Mn = min(Mn, min(Mn_LTB, Mp))
        end
    end

    # FLB (F3-1 / F3-2)
    if λ_f > λr_f
        kc = clamp(4 / sqrt(λ_w), 0.35, 0.76)
        Mn = min(Mn, 0.9 * E * kc * Sx / λ_f^2)
    elseif λ_f > λp_f
        Mn = min(Mn, Mp - (Mp - 0.7 * Fy * Sx) * ((λ_f - λp_f) / (λr_f - λp_f)))
    end

    return Mn
end

# ============================================================================
# Nominal flexural strength — Weak axis (F6)
# ============================================================================

function _compute_Mn_weak(Zy, Sy, λ_f, E, Fy)
    Mp = min(Fy * Zy, 1.6 * Fy * Sy)
    Mn = Mp

    λp = 0.38 * sqrt(E / Fy)
    λr = 1.0 * sqrt(E / Fy)

    if λ_f > λr
        Fcr = 0.69 * E / λ_f^2
        Mn = min(Mn, Fcr * Sy)
    elseif λ_f > λp
        Mn = min(Mn, Mp - (Mp - 0.7 * Fy * Sy) * ((λ_f - λp) / (λr - λp)))
    end

    return Mn
end

# ============================================================================
# Nominal shear strength — Strong axis (G2.1)
# ============================================================================

function _compute_Vn_strong(d, tw, tf, E, Fy; kv=5.34)
    Aw = d * tw
    h_w = d - 2 * tf

    limit_rolled = 2.24 * sqrt(E / Fy)
    if h_w / tw <= limit_rolled
        Cv1 = 1.0
    else
        limit_inelastic = 1.10 * sqrt(kv * E / Fy)
        if h_w / tw <= limit_inelastic
            Cv1 = 1.0
        else
            Cv1 = limit_inelastic / (h_w / tw)
        end
    end

    return 0.6 * Fy * Aw * Cv1
end

# ============================================================================
# Nominal shear strength — Weak axis (G7)
# ============================================================================

function _compute_Vn_weak(bf, tf, Fy)
    Aw_weak = 2 * bf * tf
    Cv2 = 1.0  # compact flanges for rolled shapes
    return 0.6 * Fy * Aw_weak * Cv2
end

# ============================================================================
# Section coordinates for visualization
# ============================================================================

function _compute_coords(d, bf, tw, tf)
    return [
        [-bf/2, 0],
        [ bf/2, 0],
        [ bf/2, -tf],
        [ tw/2, -tf],
        [ tw/2, -(d-tf)],
        [ bf/2, -(d-tf)],
        [ bf/2, -d],
        [-bf/2, -d],
        [-bf/2, -(d-tf)],
        [-tw/2, -(d-tf)],
        [-tw/2, -tf],
        [-bf/2, -tf],
        [-bf/2, 0]
    ]
end

# ============================================================================
# Mutators
# ============================================================================

"""
    update!(self::I_symm; d, bf, tw, tf, material, Lb)

Modify section in place and recompute all properties.
"""
function update!(self::I_symm; d::Real=self.d, bf::Real=self.bf,
                 tw::Real=self.tw, tf::Real=self.tf,
                 material::Metal=self.material, Lb::Real=self.Lb)
    setfield!(self, :d, d)
    setfield!(self, :bf, bf)
    setfield!(self, :tw, tw)
    setfield!(self, :tf, tf)
    setfield!(self, :material, material)
    setfield!(self, :Lb, Lb)
    get_properties!(self)
end

function update!(self::I_symm, vars::Vector{T}) where T <: Real
    update!(self, d=vars[1], bf=vars[2], tw=vars[3], tf=vars[4])
end

"""
    update(self::I_symm; d, bf, tw, tf, material, Lb) -> I_symm

Return a new section with modified parameters.
"""
function update(self::I_symm; d::Real=self.d, bf::Real=self.bf,
                tw::Real=self.tw, tf::Real=self.tf,
                material::Metal=self.material, Lb::Real=self.Lb)
    return I_symm(d, bf, tw, tf, material=material, Lb=Lb)
end

# ============================================================================
# Geometry extraction
# ============================================================================

"""
    get_geometry_vars(self::I_symm) -> Vector{Real}

Return `[d, bf, tw, tf]` for use in optimization and catalog lookup.
"""
get_geometry_vars(self::I_symm) = [self.d, self.bf, self.tw, self.tf]

"""
    get_coords(self::I_symm) -> Vector{Vector{Float64}}

Return 2D outline coordinates of the cross-section.
"""
get_coords(self::I_symm) = getfield(self, :flat_coords)

# ============================================================================
# LRFD design capacities — convenience functions
# ============================================================================

"""
    get_Mn(self::I_symm; Lb, Cb, axis) -> Real

Nominal flexural strength per AISC 360-16 Chapter F.
Recomputes for the given `Lb`/`Cb` without mutating the struct.
"""
function get_Mn(self::I_symm; Lb::Real=self.Lb, Cb::Real=1.0, axis::Symbol=:strong)
    E, Fy = self.material.E, self.material.Fy
    if axis == :strong
        if Lb == getfield(self, :Lb) && Cb == 1.0
            return getfield(self, :Mn)
        end
        Lp_val, Lr_val, Fcr_val = if Lb > 0
            _compute_LTB_params(self.ry, self.J, self.Sx, self.ho, self.rts,
                                self.Cw, self.Iy, Lb, E, Fy)
        else
            (0.0, 0.0, 0.0)
        end
        return _compute_Mn_strong(self.Zx, self.Sx, self.λ_f, self.λp_f, self.λr_f,
                                  self.λ_w, Lb, Lp_val, Lr_val, Fcr_val, E, Fy)
    else
        return getfield(self, :Mn_weak)
    end
end

"""
    get_ϕMn(self::I_symm; Lb, Cb, axis, ϕ) -> Real

Design flexural strength ϕMn (LRFD, ϕ_b = 0.9).
"""
get_ϕMn(self::I_symm; Lb::Real=self.Lb, Cb::Real=1.0,
        axis::Symbol=:strong, ϕ::Real=0.9) =
    ϕ * get_Mn(self; Lb=Lb, Cb=Cb, axis=axis)

"""
    get_Vn(self::I_symm; axis, kv) -> Real

Nominal shear strength per AISC 360-16 Chapter G.
"""
function get_Vn(self::I_symm; axis::Symbol=:strong, kv::Real=5.34)
    if axis == :strong
        return getfield(self, :Vn)
    else
        return getfield(self, :Vn_weak)
    end
end

"""
    get_ϕVn(self::I_symm; axis, ϕ) -> Real

Design shear strength ϕVn (LRFD). Default ϕ = 1.0 for strong-axis
rolled I-shapes per G2.1(a), ϕ = 0.9 for weak-axis per G7.
"""
function get_ϕVn(self::I_symm; axis::Symbol=:strong, ϕ::Union{Real,Nothing}=nothing)
    ϕ_use = isnothing(ϕ) ? (axis == :strong ? 1.0 : 0.9) : ϕ
    return ϕ_use * get_Vn(self; axis=axis)
end

# ============================================================================
# Compression (Chapter E)
# ============================================================================

"""
    get_compression_factors(self::I_symm) -> NamedTuple(:Qs, :Qa, :Q)

Slender element reduction factors per AISC 360-16 E7.
"""
function get_compression_factors(self::I_symm)
    E, Fy = self.material.E, self.material.Fy
    λ_f = getfield(self, :λ_f)
    λ_w = getfield(self, :λ_w)

    # Qs: Unstiffened (flanges) — rolled I-shapes
    qs_lim1 = 0.56 * sqrt(E / Fy)
    qs_lim2 = 1.03 * sqrt(E / Fy)
    if λ_f <= qs_lim1
        Qs = 1.0
    elseif λ_f < qs_lim2
        Qs = 1.415 - 0.74 * λ_f * sqrt(Fy / E)
    else
        Qs = 0.69 * E / (Fy * λ_f^2)
    end

    # Qa: Stiffened (webs) — E7 effective width
    λr_w_comp = 1.49 * sqrt(E / Fy)
    if λ_w <= λr_w_comp
        Qa = 1.0
    else
        h_w = getfield(self, :d) - 2 * getfield(self, :tf)
        c1, c2 = 0.18, 1.31
        Fel = (c2 * λr_w_comp / λ_w)^2 * Fy
        Fcr_comp = Fy
        ratio = sqrt(Fel / Fcr_comp)
        be = clamp(h_w * ratio * (1 - c1 * ratio), 0.0, h_w)
        Ae = getfield(self, :A) - (h_w - be) * getfield(self, :tw)
        Qa = max(Ae / getfield(self, :A), 0.0)
    end

    return (Qs=Qs, Qa=Qa, Q=Qs * Qa)
end

"""
    get_Fe_flexural(self::I_symm, L; K=1.0, axis=:weak) -> Real

Elastic flexural buckling stress Fe = π²E/(KL/r)² per AISC 360-16 E3-4.

`K` is the effective length factor for the given axis (Kx for strong, Ky for weak).
"""
function get_Fe_flexural(self::I_symm, L::Real; K::Real=1.0, axis::Symbol=:weak)
    r = axis == :weak ? self.ry : self.rx
    return π^2 * self.material.E / (K * L / r)^2
end

"""
    get_Fe_torsional(self::I_symm, Lz; Kz=1.0) -> Real

Elastic torsional buckling stress per AISC 360-16 E4-4.
"""
function get_Fe_torsional(self::I_symm, Lz::Real; Kz::Real=1.0)
    E, G = self.material.E, self.material.G
    return (π^2 * E * self.Cw / (Kz * Lz)^2 + G * self.J) / (self.Ix + self.Iy)
end

"""
    get_Pn(self, Lx, Ly; Kx=1.0, Ky=1.0) -> Real
    get_Pn(self, L; axis=:weak)            -> Real (legacy)

Nominal compressive strength per AISC 360-16 E3/E4/E7.

The two-length form computes Fe for both axes and takes the minimum (governing).
The single-length form is kept for backward compatibility.

# Keyword arguments
- `Kx`, `Ky` : effective length factors (strong / weak axis)
"""
function get_Pn(self::I_symm, Lx::Real, Ly::Real;
                Kx::Real=1.0, Ky::Real=1.0)
    Q  = get_compression_factors(self).Q
    Fy = self.material.Fy

    Fe_x = get_Fe_flexural(self, Lx; K=Kx, axis=:strong)
    Fe_y = get_Fe_flexural(self, Ly; K=Ky, axis=:weak)
    Fe_t = get_Fe_torsional(self, max(Lx, Ly))
    Fe   = min(Fe_x, Fe_y, Fe_t)

    ratio = Q * Fy / Fe
    Fcr = ratio <= 2.25 ?
          Q * 0.658^ratio * Fy :
          0.877 * Fe

    return Fcr * self.A
end

function get_Pn(self::I_symm, L::Real; axis::Symbol=:weak)
    Q = get_compression_factors(self).Q
    Fy = self.material.Fy

    Fe = axis == :torsional ?
         get_Fe_torsional(self, L) :
         get_Fe_flexural(self, L; axis=axis)

    ratio = Q * Fy / Fe
    Fcr = ratio <= 2.25 ?
          Q * 0.658^ratio * Fy :
          0.877 * Fe

    return Fcr * self.A
end

"""
    get_ϕPn(self, Lx, Ly; Kx=1.0, Ky=1.0, ϕ=0.90) -> Real
    get_ϕPn(self, L; axis=:weak, ϕ=0.90)            -> Real (legacy)

Design compressive strength ϕPn (LRFD, ϕ_c = 0.90).
"""
get_ϕPn(self::I_symm, Lx::Real, Ly::Real;
        Kx::Real=1.0, Ky::Real=1.0, ϕ::Real=0.90) =
    ϕ * get_Pn(self, Lx, Ly; Kx=Kx, Ky=Ky)

get_ϕPn(self::I_symm, L::Real; axis::Symbol=:weak, ϕ::Real=0.90) =
    ϕ * get_Pn(self, L; axis=axis)

# ============================================================================
# Slenderness queries
# ============================================================================

"""
    get_slenderness(self::I_symm) -> NamedTuple

Return flange/web slenderness ratios, limits, and classifications.
"""
function get_slenderness(self::I_symm)
    return (λ_f=self.λ_f, λp_f=self.λp_f, λr_f=self.λr_f,
            class_f=getfield(self, :slenderness_f),
            λ_w=self.λ_w, λp_w=self.λp_w, λr_w=self.λr_w,
            class_w=getfield(self, :slenderness_w))
end

"""
    is_compact(self::I_symm) -> Bool

Check whether the section is compact in both flange and web for flexure.
"""
function is_compact(self::I_symm)
    return getfield(self, :slenderness_f) == :compact &&
           getfield(self, :slenderness_w) == :compact
end

# Keep old name as alias
get_lambdas(self::I_symm) = (self.λ_f, self.λp_f, self.λr_f,
                              self.λ_w, self.λp_w, self.λr_w)

get_Ls(self::I_symm) = (getfield(self, :Lp), getfield(self, :Lr), getfield(self, :F_cr))

# ============================================================================
# Constraint vector for optimization
# ============================================================================

"""
    get_constraint_vector(d, bf, tw, tf, Mu, Vu; ϕ_b=0.9, ϕ_v=0.9) -> Vector

Return a vector of inequality constraints (each ≤ 0 when feasible) for use in
NLopt or sequential search.

Constraints: [-A, 2tf-d, Aw/Afc-10, h_w/tw-260, 0.1-Iyc/Iy, Iyc/Iy-0.9, moment util, shear util]
"""
function get_constraint_vector(d::Real, bf::Real, tw::Real, tf::Real,
                               Mu_section::Real, Vu_section::Real;
                               ϕ_b::Real=0.9, ϕ_v::Real=0.9)
    section = I_symm(d, bf, tw, tf)

    A     = -section.A
    Δ_height = 2 * tf - d

    h_w   = d - 2 * tf
    A_w   = tw * h_w
    A_fc  = tf * bf
    ratio_w_to_f  = A_w / A_fc - 10
    ratio_d_to_tw = h_w / tw - 260

    flange_ratio_min = 0.1 - section.I_yc / section.Iy
    flange_ratio_max = section.I_yc / section.Iy - 0.9

    util_moment = abs(Mu_section) / abs(ϕ_b * section.Mn) - 1.0
    util_shear  = abs(Vu_section) / abs(ϕ_v * section.Vn) - 1.0

    return [A, Δ_height, ratio_w_to_f, ratio_d_to_tw,
            flange_ratio_min, flange_ratio_max,
            util_moment, util_shear]
end
