"""
Composite beam stiffness for serviceability deflection checks.

All units are inches and ksi (consistent with VariableBeamOptimizer conventions).
Implements AISC I3.1a effective width and elastic transformed-section moment of inertia
for full-composite solid slab on steel I-beam.

Self-contained: inlines the Ix calculation so this file has no external dependencies.
"""

"""
    _Ix_I_symm(h, w, tw, tf)

Strong-axis moment of inertia for a symmetric I/H section (inlined from
VariableBeamOptimizer to avoid include-order dependency).
"""
function _Ix_I_symm(h::Real, w::Real, tw::Real, tf::Real)
    h_w = h - 2 * tf
    A_b = w * tf
    A_w = h_w * tw
    A_t = w * tf
    A = A_b + A_w + A_t

    y_b = tf / 2
    y_w = tf + h_w / 2
    y_t = h - tf / 2
    y_c = (A_b * y_b + A_w * y_w + A_t * y_t) / A

    I_b = w * tf^3 / 12 + A_b * (y_b - y_c)^2
    I_w = tw * h_w^3 / 12 + A_w * (y_w - y_c)^2
    I_t = w * tf^3 / 12 + A_t * (y_t - y_c)^2

    return I_b + I_w + I_t
end

"""
    get_b_eff(L_beam, trib_width_total; is_perimeter=false)

Effective concrete flange width per AISC 360, Section I3.1a.

Works directly from tributary polygon geometry — no "beam spacing" needed.

**Interior beams** (slab on both sides): `trib_width_total` is the sum of
the one-sided tributary widths from the two adjacent cycles. Each side
contributes at most `L/8`, so `b_eff = min(L/4, trib_width_total)`.

**Perimeter beams** (slab on one side only): only one cycle contributes
tributary width. The slab-side half-width is capped at `L/8`, and there
is no slab on the free edge, so `b_eff = min(L/8, trib_width_total)`.

# Arguments
- `L_beam`: beam span (in)
- `trib_width_total`: total perpendicular tributary width at this station (in),
  already aggregated across both sides for interior beams
- `is_perimeter`: true when the beam is on the slab perimeter (slab one side only)
"""
function get_b_eff(L_beam::Real, trib_width_total::Real; is_perimeter::Bool=false)
    if is_perimeter
        return min(L_beam / 8, trib_width_total)
    else
        return min(L_beam / 4, trib_width_total)
    end
end

"""
    get_I_composite(h, w, tw, tf, t_slab, b_eff, E_s, E_c)

Full-composite elastic (transformed) moment of inertia for a solid concrete
slab on a symmetric steel I-section. Used for serviceability deflection checks.

The concrete slab is transformed to an equivalent steel width via the modular
ratio n = E_s / E_c, then the neutral axis and parallel-axis theorem are applied.

# Arguments
- `h, w, tw, tf`: steel I-section dimensions (in) — depth, flange width, web thickness, flange thickness
- `t_slab`: concrete slab thickness (in)
- `b_eff`: effective slab width (in), from `get_b_eff`
- `E_s`: steel elastic modulus (ksi)
- `E_c`: concrete elastic modulus (ksi)

# Returns
- Transformed composite moment of inertia I_tr (in⁴)
"""
function get_I_composite(h::Real, w::Real, tw::Real, tf::Real,
                         t_slab::Real, b_eff::Real, E_s::Real, E_c::Real)

    n = E_s / E_c

    h_w = h - 2 * tf
    A_steel = 2 * w * tf + h_w * tw

    b_tr = b_eff / n
    A_conc_tr = b_tr * t_slab

    y_steel = h / 2
    y_conc = h + t_slab / 2

    A_total = A_steel + A_conc_tr
    y_bar = (A_steel * y_steel + A_conc_tr * y_conc) / A_total

    I_steel = _Ix_I_symm(h, w, tw, tf)

    I_conc_own = b_tr * t_slab^3 / 12

    I_tr = I_steel + A_steel * (y_steel - y_bar)^2 +
           I_conc_own + A_conc_tr * (y_conc - y_bar)^2

    return I_tr
end

"""
    composite_stiffness_ratio(h, w, tw, tf, t_slab, b_eff, E_s, E_c)

Ratio of composite I to bare steel I. Useful for understanding the stiffness
gain from composite action.
"""
function composite_stiffness_ratio(h::Real, w::Real, tw::Real, tf::Real,
                                   t_slab::Real, b_eff::Real, E_s::Real, E_c::Real)
    I_bare = _Ix_I_symm(h, w, tw, tf)
    I_comp = get_I_composite(h, w, tw, tf, t_slab, b_eff, E_s, E_c)
    return I_comp / I_bare
end

"""
    _aggregate_strips(load_positions, trib_widths; atol=1e-6)

Merge load strips that share the same beam station (within `atol`). For
interior beams that receive tributary load from two adjacent cycles, two
strips at the same station carry one-sided widths; summing them gives the
total slab width available for composite action at that station.

Returns `(agg_positions, agg_widths)` with duplicates merged.
"""
function _aggregate_strips(load_positions::Vector{<:Real},
                           trib_widths::Vector{<:Real};
                           atol::Real=1e-6)
    n = length(load_positions)
    n == 0 && return (Float64[], Float64[])

    sp = sortperm(load_positions)
    agg_pos = Float64[load_positions[sp[1]]]
    agg_w   = Float64[trib_widths[sp[1]]]

    for k in 2:n
        if abs(load_positions[sp[k]] - agg_pos[end]) <= atol
            agg_w[end] += trib_widths[sp[k]]
        else
            push!(agg_pos, load_positions[sp[k]])
            push!(agg_w, trib_widths[sp[k]])
        end
    end
    return (agg_pos, agg_w)
end

"""
    get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                              L_beam, load_positions, trib_widths;
                              is_perimeter=false)

Effective composite moment of inertia for a beam with spatially varying
tributary width, computed via virtual-work weighting so that regions of
high curvature (near midspan) dominate the effective stiffness.

For a simply-supported beam the virtual moment shape from a unit midspan
load is M̄(x) = x(L−x)/(2L).  The effective I is defined by

    1/I_eff = [∫ M̄(x)²/I(x) dx] / [∫ M̄(x)² dx]

where I(x) = I_composite at the local tributary width.

**Strip aggregation**: load points at the same beam station (from two
adjacent cycles) are merged by summing their one-sided tributary widths.

**Perimeter beams**: when `is_perimeter=true`, the slab extends only
on one side, so the AISC single-side cap applies: `b_eff = min(L/8, w)`
instead of the interior `min(L/4, w)`.

Falls back to bare steel Ix when no load points are supplied.

# Arguments
- `h, w, tw, tf`: steel section dimensions (in)
- `t_slab`: slab thickness (in)
- `E_s, E_c`: elastic moduli (ksi)
- `L_beam`: beam span (in)
- `load_positions`: normalized positions along beam (0–1) of each strip
- `trib_widths`: perpendicular tributary width (in) at each strip
- `is_perimeter`: true for edge beams (slab on one side only)
"""
function get_I_composite_effective(h::Real, w::Real, tw::Real, tf::Real,
                                   t_slab::Real, E_s::Real, E_c::Real,
                                   L_beam::Real,
                                   load_positions::Vector{<:Real},
                                   trib_widths::Vector{<:Real};
                                   is_perimeter::Bool=false)

    n_raw = length(load_positions)
    n_raw == 0 && return _Ix_I_symm(h, w, tw, tf)

    agg_pos, agg_w = _aggregate_strips(load_positions, trib_widths)
    n_pts = length(agg_pos)

    x_abs = agg_pos .* L_beam

    M_bar = x_abs .* (L_beam .- x_abs)

    sp = sortperm(x_abs)
    x_sorted = x_abs[sp]

    dx = zeros(n_pts)
    if n_pts == 1
        dx[sp[1]] = L_beam
    else
        for (idx, si) in enumerate(sp)
            if idx == 1
                dx[si] = (x_sorted[2] - x_sorted[1]) / 2 + x_sorted[1]
            elseif idx == n_pts
                dx[si] = (L_beam - x_sorted[end]) + (x_sorted[end] - x_sorted[end-1]) / 2
            else
                dx[si] = (x_sorted[idx+1] - x_sorted[idx-1]) / 2
            end
        end
    end

    I_local = Vector{Float64}(undef, n_pts)
    for k in 1:n_pts
        b_eff_k = get_b_eff(L_beam, agg_w[k]; is_perimeter=is_perimeter)
        I_local[k] = get_I_composite(h, w, tw, tf, t_slab, b_eff_k, E_s, E_c)
    end

    num = sum(M_bar[k]^2 * dx[k] for k in 1:n_pts)
    den = sum(M_bar[k]^2 / I_local[k] * dx[k] for k in 1:n_pts)

    (num <= 0 || den <= 0) && return _Ix_I_symm(h, w, tw, tf)
    return num / den
end
