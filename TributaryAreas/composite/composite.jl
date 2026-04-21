"""
Composite beam stiffness for serviceability deflection checks.

All units are inches and ksi (consistent with VariableBeamOptimizer conventions).
Implements AISC I3.1a effective width and elastic transformed-section moment of inertia
for full-composite solid slab on steel I-beam.

Self-contained: inlines the Ix calculation so this file has no external
dependencies. If Zygote is loaded in the including scope (the normal case
when running from `_scripts.jl`), `get_I_composite_effective` hides its
non-differentiable geometry preprocessing from AD; otherwise the helper
below becomes a no-op and the function still runs normally.
"""

"""
    _ad_ignore(f)

Run `f()` with Zygote AD hidden if Zygote is available in the including
module, otherwise just call `f()`. Keeps `composite.jl` usable both from
the full optimization stack (where Zygote is required) and from
standalone diagnostics like `test_composite_effective.jl`.
"""
@inline function _ad_ignore(f)
    if @isdefined(Zygote)
        return Zygote.ignore(f)
    else
        return f()
    end
end

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
the per-side tributary widths from the two adjacent cycles, *each already
capped at `L/8`* by the caller (see `get_I_composite_effective`). The
total is therefore bounded by `L/4`, matching AISC I3.1a(a).

**Perimeter beams** (slab on one side only): only one cycle contributes
tributary width. The slab-side half-width is capped at `L/8`, and there
is no slab on the free edge, so `b_eff = min(L/8, trib_width_total)`.

Applying the per-side L/8 cap *before* summing matters for irregular
geometries where one side's trib is much wider than the other: summing
raw widths and capping the sum at L/4 lets a wide side compensate for a
narrow side, overestimating the effective flange. See the comments in
`get_I_composite_effective` for the pre-aggregation handling.

# Arguments
- `L_beam`: beam span (in)
- `trib_width_total`: total perpendicular tributary width at this station (in),
  already per-side capped and aggregated across both sides for interior beams
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
tributary width, computed as a moment-squared-weighted harmonic mean of
the local composite `I(x)` so that regions of high curvature (near
midspan) dominate the effective stiffness.

For a simply-supported beam under a uniformly distributed load the
applied moment shape is `M(x) ∝ x(L−x)`. The effective I is defined by

    1/I_eff = [∫ M(x)²/I(x) dx] / [∫ M(x)² dx]

which is the strain-energy-weighted reciprocal stiffness (equivalent to
a midspan-deflection virtual-work integral up to the difference between
the applied-moment and unit-load-moment shapes — that difference is
small and conservative for the typical trib profiles seen here).

**Per-side L/8 cap**: each strip carries a *one-sided* tributary width
from a single cycle. The AISC I3.1a cap is L/8 *per side*, so each strip
width is capped at L/8 **before** the two-side aggregation. Summing raw
widths and only capping the total at L/4 would let a wide side
compensate for a narrow side.

**Strip aggregation**: after per-side capping, load points at the same
beam station (from two adjacent cycles bounding an interior beam) are
merged by summing their capped one-sided widths.

**Perimeter beams**: when `is_perimeter=true`, the slab extends only
on one side, so only one cycle contributes a strip per station. The
L/8 per-side cap is then the controlling limit (`get_b_eff` enforces
`min(L/8, w)`).

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

    # -------------------------------------------------------------------------
    # Load-geometry preprocessing.
    #
    # `load_positions`, `trib_widths`, and `L_beam` are constants from the
    # caller's perspective — none depend on the steel section dimensions
    # (h, w, tw, tf). We therefore hide the whole preprocessing block from
    # Zygote so the `push!`/`setindex!`/`sortperm` patterns in
    # `_aggregate_strips` and the neighbour-difference build of `dx` don't
    # break AD when this function is called from the NLP deflection
    # constraint. Results returned as `Vector{Float64}` triples.
    #
    # Per-side L/8 cap (AISC I3.1a(a)): each strip represents one side's
    # tributary; capping here enforces the per-side limit, so the total
    # across two cycles is bounded by L/4 for interior beams. Perimeter
    # beams carry only one side per station, so the single-side L/8 cap
    # is controlling.
    # -------------------------------------------------------------------------
    agg_w_sorted, M_shape_sq_sorted, dx_sorted = _ad_ignore() do
        side_cap             = L_beam / 8
        trib_widths_capped   = min.(trib_widths, side_cap)
        agg_pos, agg_w_local = _aggregate_strips(load_positions, trib_widths_capped)
        n_pts_local          = length(agg_pos)

        x_abs       = agg_pos .* L_beam
        sp          = sortperm(x_abs)
        x_sorted    = x_abs[sp]
        agg_w_s     = agg_w_local[sp]

        # Applied-moment shape under a UDL on a simply-supported beam,
        # pre-squared for the energy weighting used in num/den below.
        M_sq_s = (x_sorted .* (L_beam .- x_sorted)) .^ 2

        dx_s = if n_pts_local == 1
            Float64[L_beam]
        else
            [
                (k == 1)            ? (x_sorted[2] - x_sorted[1]) / 2 + x_sorted[1] :
                (k == n_pts_local)  ? (L_beam - x_sorted[end]) + (x_sorted[end] - x_sorted[end-1]) / 2 :
                                      (x_sorted[k+1] - x_sorted[k-1]) / 2
                for k in 1:n_pts_local
            ]
        end

        return agg_w_s, M_sq_s, dx_s
    end

    n_pts = length(agg_w_sorted)

    # -------------------------------------------------------------------------
    # Gradient-tracked: I_local depends on the section dimensions through
    # `get_I_composite`. Built as a comprehension (no setindex!) so Zygote
    # can propagate the pullback back to (h, w, tw, tf).
    # -------------------------------------------------------------------------
    I_local = [get_I_composite(h, w, tw, tf, t_slab,
                               get_b_eff(L_beam, agg_w_sorted[k]; is_perimeter=is_perimeter),
                               E_s, E_c)
               for k in 1:n_pts]

    num = sum(M_shape_sq_sorted[k] * dx_sorted[k] for k in 1:n_pts)
    den = sum(M_shape_sq_sorted[k] / I_local[k] * dx_sorted[k] for k in 1:n_pts)

    (num <= 0 || den <= 0) && return _Ix_I_symm(h, w, tw, tf)
    return num / den
end
