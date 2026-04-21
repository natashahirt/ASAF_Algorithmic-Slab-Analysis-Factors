"""
Head-to-head comparison of the Nov 2024 bare-steel two-stage NLP sizer
(`size_bare_nlp!` from `VariableBeamOptimizer/sizer_nlp_v1/`) against the
current composite-aware sizer (`optimal_beamsizer` with
`composite_action=true`).

Scientific question
-------------------
The paper's BaU baseline (kg CO₂eq / m² per topology) was produced by the
Nov 2024 sizer, which applied a strict **bare-steel L/360** deflection
constraint — no composite stiffness credit. When evaluated against the
simplified single-limit composite pipeline (total unfactored load on the
AISC §I3.1a transformed section, checked to L/360 in one shot), are
those Nov 2024 sections:

  1. **conservative** on mass (`mass_nlp ≥ mass_comp`), and
  2. **feasible** under composite strength (M_u ≤ ϕM_n) and single-limit
     composite deflection (L/360 on transformed Ix)?

A "yes" on both validates the paper's embodied-carbon figures without
re-running the 2024 sweeps.

Per-bay protocol
----------------
For each topology JSON in the scope list (BaU by default, expandable to
any `name ∈ Geometries/topology/`):

  1. **Size A (`nlp`)** — call `size_bare_nlp!` with the Nov 2024 regime
     (`composite_action=false`, `deflection_reduction_factor=1.0`,
     `minimum_continuous=true`, W4X13 floor, discrete-then-continuous
     warm start). Record `mass_nlp = Σ Aᵢ · Lᵢ · ρ_steel`. This is the
     primary cost (~4 min per 174-beam bay).

  2. **Evaluate A under composite FE** — build a fresh
     `SlabSizingParams` with `composite_action=true`, copy the Size-A
     minimizers into it, and call `postprocess_slab`. Harvest
     `strength_ok`, `n_L360_fail`, `max_util_M`, etc. This is the
     modern-composite verdict on the Nov 2024 sections and is the
     main scientific question of the comparison (~10 s/bay).

  3. **Size B (`comp`)** — only when `COMPARE_INCLUDE_COMPOSITE=1`.
     Calls `optimal_beamsizer` with `composite_action=true`,
     `beam_sizer=:continuous`, and the Size-A minimizers as
     `initial_vars`. Single-pass, single-limit regime (see
     `optimal_beamsizer`). Reports `mass_comp` and
     `mass_ratio = mass_nlp / mass_comp`.

Emits `per_bay.csv` with one row per topology and a concise summary
table to stdout.

Usage
-----
    # Scope = :bau  →  one-bay sanity check (r6c4), ~5 min
    julia --project=. SlabDesignFactors/scripts/compare_nlp_vs_composite.jl

    # Scope = :medium  →  BaU + 10 random bays, ~50 min (NLP + eval only)
    COMPARE_SCOPE=medium julia --project=. SlabDesignFactors/scripts/compare_nlp_vs_composite.jl

    # Scope = :all     →  every topology under Geometries/topology/. Long-running.
    COMPARE_SCOPE=all julia --project=. SlabDesignFactors/scripts/compare_nlp_vs_composite.jl

    # Add composite re-sizing (reconciliation loop; tightens against
    # postprocess until δ/Δ_lim ≤ 1 + tol)
    COMPARE_INCLUDE_COMPOSITE=1 julia --project=. SlabDesignFactors/scripts/compare_nlp_vs_composite.jl

    # Partial-composite-action knockdown (1.0 = full composite per AISC
    # I3.1a; 0.85 = typical partial shear connection / slip / creep)
    COMPARE_PCF=0.85 COMPARE_INCLUDE_COMPOSITE=1 julia --project=. ...

    # Max reconciliation iterations (2–3 is enough; 0 disables)
    COMPARE_RECONCILE_ITERS=3 COMPARE_INCLUDE_COMPOSITE=1 julia --project=. ...
"""

include("_scripts.jl")

using CSV, DataFrames, JSON, Random, Statistics

# ── paths / configuration ────────────────────────────────────────────────────
const GEOM_DIR    = joinpath(@__DIR__, "..", "..", "Geometries", "topology")
const RESULTS_DIR = joinpath(@__DIR__, "..", "results", "compare_nlp_vs_composite")
const LEGACY_CSV  = joinpath(@__DIR__, "..", "results",
                             "remote_results_yesdeflection_yesslabmin",
                             "topology.csv")

const SCOPE             = Symbol(get(ENV, "COMPARE_SCOPE", "bau"))  # :bau | :medium | :all
const INCLUDE_COMPOSITE = get(ENV, "COMPARE_INCLUDE_COMPOSITE", "1") == "1"
const MAX_DEPTH_IN      = 25.0
const SEED              = 42

# Partial-composite-action knockdown on Ix_composite, shared between the
# sizer constraint and the postprocess FE so the two paths stay consistent.
# 1.0 = full AISC I3.1a transformed section; 0.85 is a typical value for
# partial shear connection, slip, and long-term concrete effects.
const PARTIAL_COMPOSITE_FACTOR = parse(Float64, get(ENV, "COMPARE_PCF", "1.0"))

# Reconciliation loop controls: if the postprocess FE reports δ/Δ_lim > 1,
# multiply the sizer's `serviceability_tighten` by the measured ratio,
# warm-start from the last minimizers, and re-size. 0 disables the loop.
const RECONCILE_MAX_ITER = parse(Int,     get(ENV, "COMPARE_RECONCILE_ITERS", "2"))
const RECONCILE_TOL      = parse(Float64, get(ENV, "COMPARE_RECONCILE_TOL",   "0.05"))

mkpath(RESULTS_DIR)

# ── bay scope definitions ────────────────────────────────────────────────────
"""
    _bau_row() -> NamedTuple

BaU design point used throughout the paper. Mirrors the reproduction test so
the NLP mass is directly comparable to its Nov 2024 CSV entry.
"""
_bau_row() = (
    name        = "r6c4",
    slab_type   = :isotropic,
    slab_sizer  = :uniform,
    collinear   = false,
    vector_1d   = [0.0, 0.0],
    max_depth   = MAX_DEPTH_IN,
)

"""
    _random_bays(n; seed) -> Vector{NamedTuple}

Draw `n` distinct topology names (BaU excluded) from the Geometries/topology/
folder, each with the BaU design point (isotropic / uniform / non-collinear).
"""
function _random_bays(n::Integer; seed::Integer=SEED)
    all_names = [splitext(basename(p))[1] for p in readdir(GEOM_DIR; join=true)
                 if endswith(p, ".json")]
    filter!(!=("r6c4"), all_names)
    rng = MersenneTwister(seed)
    n_draw = min(n, length(all_names))
    picks  = all_names[randperm(rng, length(all_names))[1:n_draw]]
    return [(; name,
              slab_type = :isotropic,
              slab_sizer = :uniform,
              collinear  = false,
              vector_1d  = [0.0, 0.0],
              max_depth  = MAX_DEPTH_IN) for name in picks]
end

"""
    _all_bays() -> Vector{NamedTuple}

Every topology JSON in `Geometries/topology/`, each with the BaU design
point. Use sparingly — running the full NLP + reconciled composite sweep
across all topologies is hours of wall time.
"""
function _all_bays()
    names = [splitext(basename(p))[1] for p in readdir(GEOM_DIR; join=true)
             if endswith(p, ".json")]
    # BaU first so the most familiar result appears early in the CSV.
    bau = "r6c4"
    if bau in names
        filter!(!=(bau), names)
        pushfirst!(names, bau)
    end
    return [(; name,
              slab_type = :isotropic,
              slab_sizer = :uniform,
              collinear  = false,
              vector_1d  = [0.0, 0.0],
              max_depth  = MAX_DEPTH_IN) for name in names]
end

function _scope_bays(scope::Symbol)
    scope === :bau    && return [_bau_row()]
    scope === :medium && return vcat([_bau_row()], _random_bays(10))
    scope === :all    && return _all_bays()
    error("Unknown scope $scope — use :bau, :medium, or :all.")
end

# ── per-bay helpers ──────────────────────────────────────────────────────────
"""
    _load_geom(name) -> Dict

Load a topology JSON into a `Dict`, matching the sweep driver's doubled-parse
idiom for escape-sequence handling.
"""
function _load_geom(name::AbstractString)
    path = joinpath(GEOM_DIR, name * ".json")
    raw  = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""),
                                  dicttype=Dict))
    return raw isa Dict ? raw : Dict(pairs(raw))
end

"""
    _make_slab_params(geom_dict, bay) -> SlabAnalysisParams

Build the `SlabAnalysisParams` shared by both sizing regimes.
"""
function _make_slab_params(geom_dict, bay)
    geom, _ = Base.invokelatest(generate_from_json, geom_dict;
                                plot=false, drawn=false)
    return SlabAnalysisParams(
        geom,
        slab_name     = bay.name,
        slab_type     = bay.slab_type,
        vector_1d     = bay.vector_1d,
        slab_sizer    = bay.slab_sizer,
        spacing       = 0.1,
        plot_analysis = false,
        fix_param     = true,
        slab_units    = :m,
    )
end

"""
    _make_sizing_params(; composite, beam_sizer) -> SlabSizingParams

Nov 2024 load set (50 psf live, 15 psf SDL, 1.6/1.2 factors, L/360,
`minimum_continuous=true`). `composite` toggles the modern composite-aware
pipeline; `beam_sizer=:continuous` forces the NLP path.
"""
function _make_sizing_params(; composite::Bool, beam_sizer::Symbol,
                               collinear::Bool,
                               max_depth::Real=MAX_DEPTH_IN,
                               reconcile_iters::Integer=0)
    return SlabSizingParams(
        live_load              = psf_to_ksi(50),
        superimposed_dead_load = psf_to_ksi(15),
        slab_dead_load         = 0.0,
        live_factor            = 1.6,
        dead_factor            = 1.2,
        beam_sizer             = beam_sizer,
        nlp_solver             = :MMA,
        max_depth              = max_depth,
        beam_units             = :in,
        serviceability_lim     = 360,
        collinear              = collinear,
        minimum_continuous     = true,
        n_max_sections         = 0,
        composite_action       = composite,
        # Partial-composite knockdown only meaningful for composite mode;
        # bare-steel Nov 2024 path keeps the default (pcf is ignored
        # there because the sizer doesn't call `get_I_composite_effective`).
        partial_composite_factor = composite ? PARTIAL_COMPOSITE_FACTOR : 1.0,
        # Sizer-postprocess reconciliation only applies to composite mode.
        reconcile_max_iter = composite ? reconcile_iters : 0,
        reconcile_tol      = RECONCILE_TOL,
    )
end

"""
    _beam_mass_kg(sizing_params) -> Float64

Sum `A · L` over the sized beam minimizers (in-system) and convert to steel
mass in kg (same idiom as `replot_composite_mass_scatter.jl`).
"""
function _beam_mass_kg(sizing_params)
    isempty(sizing_params.minimizers) && return 0.0
    Ls = [be.length for be in sizing_params.model.elements[:beam]]
    return sum(I_symm(m...).A * Ls[i]
               for (i, m) in enumerate(sizing_params.minimizers)
              ) * convert_to_m[:in]^3 * ρ_STEEL
end

"""
    _size_nlp(bay) -> (mass_nlp, results_nlp, sizing_params_nlp)

Nov 2024 bare-steel two-stage NLP sizing. Returns the reproduced mass, the
`postprocess_slab` results (with bare-steel deflection reporting — note
this is the self-check, NOT the composite verdict), and the sizer params
for downstream composite re-evaluation.
"""
function _size_nlp(bay)
    geom_dict = _load_geom(bay.name)
    slab_params   = _make_slab_params(geom_dict, bay)
    sizing_params = _make_sizing_params(composite=false, beam_sizer=:continuous,
                                        collinear=bay.collinear,
                                        max_depth=bay.max_depth)

    slab_params                 = analyze_slab(slab_params)
    slab_params, sizing_params  = size_bare_nlp!(slab_params, sizing_params)

    results = postprocess_slab(slab_params, sizing_params;
                               check_collinear=bay.collinear)
    return _beam_mass_kg(sizing_params), results, slab_params, sizing_params
end

"""
    _size_composite(bay, nlp_minimizers) -> (mass_comp, results, n_iters, max_ratio)

Modern composite-aware **continuous** sizing under the simplified
single-limit regime: one constraint per beam,
`δ_composite ≤ L / serviceability_lim`, evaluated against the AISC
§I3.1a transformed-section `I_x` (optionally with a partial-composite
knockdown from `COMPARE_PCF`). Warm-starts from the Nov 2024 NLP
minimizers.

When `COMPARE_RECONCILE_ITERS > 0`, wraps `optimal_beamsizer` in a
postprocess-driven tightening loop (`optimal_beamsizer_reconcile`) so
the sections the sizer returns also satisfy the global-FE L/360 check
run by `postprocess_slab` — not just the per-beam sizing constraint.

Returns `(mass, results, n_iters, max_ratio)` where `results === nothing`
if the sizer produced no feasible sections.
"""
function _size_composite(bay, nlp_minimizers::Vector)
    geom_dict = _load_geom(bay.name)
    slab_params   = _make_slab_params(geom_dict, bay)
    sizing_params = _make_sizing_params(composite=true, beam_sizer=:continuous,
                                        collinear=bay.collinear,
                                        max_depth=bay.max_depth,
                                        reconcile_iters=RECONCILE_MAX_ITER)

    slab_params = analyze_slab(slab_params)

    slab_params, sizing_params, results, max_ratio, n_iters =
        optimal_beamsizer_reconcile(slab_params, sizing_params;
                                    initial_vars=copy(nlp_minimizers),
                                    verbose=true)

    if isempty(sizing_params.minimizers) || results === nothing
        return 0.0, nothing, slab_params, sizing_params, n_iters, max_ratio
    end

    # Final postprocess was already performed inside the reconciliation
    # loop. Re-run only if the caller's collinear flag differs from the
    # default used internally (currently identical, so cheap no-op).
    if bay.collinear != sizing_params.collinear
        results = postprocess_slab(slab_params, sizing_params;
                                   check_collinear=bay.collinear)
    end

    return _beam_mass_kg(sizing_params), results, slab_params, sizing_params,
           n_iters, max_ratio
end

"""
    _evaluate_nlp_under_composite(bay, nlp_minimizers, nlp_ids) -> SlabOptimResults

Build a fresh `SlabSizingParams` in composite mode, inject the Nov 2024
NLP-sized sections, and run `postprocess_slab` so the staged composite
checks (L/360 live, L/240 total, ϕM_n) fire against those sections.

This answers "do the Nov 2024 sections satisfy the modern composite
verdict?" without any re-sizing.
"""
function _evaluate_nlp_under_composite(bay, nlp_minimizers::Vector,
                                       nlp_ids::Vector{<:AbstractString})
    geom_dict = _load_geom(bay.name)
    slab_params   = _make_slab_params(geom_dict, bay)
    sizing_params = _make_sizing_params(composite=true, beam_sizer=:continuous,
                                        collinear=bay.collinear,
                                        max_depth=bay.max_depth)

    slab_params = analyze_slab(slab_params)

    # Re-scale the model and build the load dictionary the same way
    # `optimal_beamsizer` does, so `postprocess_slab` sees a consistent
    # (scaled model, load_df, load_dictionary) triple.
    sizing_params.model = slab_params.model
    conv = convert_to_m[slab_params.slab_units] * 1 / convert_to_m[sizing_params.beam_units]
    sizing_params.area  = slab_params.area * conv^2
    if sizing_params.max_assembly_depth
        slab_depth                   = maximum(slab_params.slab_depths) * conv
        sizing_params.max_beam_depth = sizing_params.max_depth - slab_depth
    end
    sizing_params.max_bay_span = maximum(slab_params.max_spans) * conv
    sizing_params.model        = get_scaled_model(slab_params, sizing_params, conv)
    sizing_params.load_dictionary = get_load_dictionary_by_id(sizing_params.model)

    # Inject the Nov 2024 NLP sections verbatim. `check_collinear=false` so
    # the diagnostics track the full beam population rather than collinear
    # representatives (the Nov 2024 CSVs are the non-collinear baseline).
    sizing_params.minimizers = deepcopy(nlp_minimizers)
    sizing_params.ids        = deepcopy(nlp_ids)
    sizing_params.minimums   = zeros(Float64, length(nlp_minimizers))

    return postprocess_slab(slab_params, sizing_params;
                            check_collinear=bay.collinear)
end

"""
    _legacy_mass(bay) -> Union{Float64,Missing}

Look up the published Nov 2024 total mass (kg) for this bay. Only the BaU
row is covered by `remote_results_yesdeflection_yesslabmin/topology.csv`
at the non-collinear / continuous / max_depth=25 configuration the bay
dict targets.
"""
function _legacy_mass(bay)
    isfile(LEGACY_CSV) || return missing
    df = CSV.read(LEGACY_CSV, DataFrame)
    mask = (df.name           .== bay.name) .&
           (string.(df.slab_type)  .== String(bay.slab_type)) .&
           (string.(df.slab_sizer) .== String(bay.slab_sizer)) .&
           (string.(df.beam_sizer) .== "continuous") .&
           (df.collinear      .== bay.collinear) .&
           (df.vector_1d_x    .≈ bay.vector_1d[1]) .&
           (df.vector_1d_y    .≈ bay.vector_1d[2]) .&
           (df.max_depth      .≈ bay.max_depth)
    rows = df[mask, :]
    nrow(rows) == 0 && return missing
    r = rows[1, :]
    return r.steel_norm * r.area
end

# ── per-beam record assembly ─────────────────────────────────────────────────
"""
    _per_beam_rows(bay, sl_nlp, sp_nlp, sl_cmp, sp_cmp) -> Vector{NamedTuple}

Build one row per beam with aligned Nov 2024 NLP and composite-sizer
geometries. Rows are keyed on `(bay_name, beam_idx)` so they can be
streamed into a flat CSV and filtered/joined downstream. Returns an
empty vector if either regime produced no feasible sections.

Columns:
  * `name`, `beam_idx`, `L_in`, `is_perimeter`
  * `A_nlp_in2` / `A_cmp_in2`      - section areas
  * `m_nlp_kg` / `m_cmp_kg`        - per-beam mass (A·L·ρ_steel)
  * `h_nlp_in` ... `tf_cmp_in`     - symmetric-I geometry (h, w, tw, tf)
"""
function _per_beam_rows(bay, sl_nlp, sp_nlp, sl_cmp, sp_cmp)
    rows = NamedTuple[]
    (isempty(sp_nlp.minimizers) || isempty(sp_cmp.minimizers)) && return rows

    # Beam-level perimeter truth lives on the slab params; the sizer only
    # populates `sizing_params.i_perimeter` in composite mode, so we
    # defer to the authoritative `slab_params.i_perimeter` for both
    # regimes.
    perim_set   = Set(sl_nlp.i_perimeter)
    beam_els    = sp_nlp.model.elements[:beam]
    n           = min(length(sp_nlp.minimizers), length(sp_cmp.minimizers),
                      length(beam_els))
    mass_factor = convert_to_m[:in]^3 * ρ_STEEL  # in³ · (kg/in³) → kg

    for i in 1:n
        h_n, w_n, tw_n, tf_n = sp_nlp.minimizers[i]
        h_c, w_c, tw_c, tf_c = sp_cmp.minimizers[i]
        sec_n = I_symm(h_n, w_n, tw_n, tf_n)
        sec_c = I_symm(h_c, w_c, tw_c, tf_c)
        L_in  = beam_els[i].length

        push!(rows, (
            name         = bay.name,
            scope        = String(SCOPE),
            beam_idx     = i,
            L_in         = L_in,
            is_perimeter = i in perim_set,
            h_nlp_in     = h_n, w_nlp_in  = w_n,
            tw_nlp_in    = tw_n, tf_nlp_in = tf_n,
            A_nlp_in2    = sec_n.A,
            m_nlp_kg     = sec_n.A * L_in * mass_factor,
            h_cmp_in     = h_c, w_cmp_in  = w_c,
            tw_cmp_in    = tw_c, tf_cmp_in = tf_c,
            A_cmp_in2    = sec_c.A,
            m_cmp_kg     = sec_c.A * L_in * mass_factor,
        ))
    end
    return rows
end

# ── per-bay runner ───────────────────────────────────────────────────────────
"""
    _compare_one(bay) -> (bay_row::NamedTuple, beam_rows::Vector{NamedTuple})

Size one bay twice, evaluate the NLP sections under composite FE, and
return both the per-bay summary and a per-beam dataframe slice keyed
on `(bay_name, beam_idx)`. `beam_rows` is empty when the composite
sizer was skipped (`COMPARE_INCLUDE_COMPOSITE=0`) or infeasible.
"""
function _compare_one(bay)
    println("\n════════════════════════════════════════════════════════════")
    println("  $(bay.name)  [scope=$SCOPE]")
    println("════════════════════════════════════════════════════════════")

    t0 = time()
    mass_nlp, res_nlp, sl_nlp, sp_nlp = _size_nlp(bay)
    t_nlp = time() - t0
    println("  NLP mass        = $(round(mass_nlp, digits=1)) kg  ($(round(t_nlp, digits=1))s)")

    legacy_mass = _legacy_mass(bay)
    if !ismissing(legacy_mass)
        println("  Legacy (2024)   = $(round(legacy_mass, digits=1)) kg")
    end

    # Evaluate the NLP sections under modern composite FE — the primary
    # scientific check. Cheap (one postprocess pass, ~10 s).
    eval_res = _evaluate_nlp_under_composite(bay,
                                             sp_nlp.minimizers, sp_nlp.ids)

    # Optional composite re-sizing (reconciled, warm-started from the NLP
    # minimizers). Off unless COMPARE_INCLUDE_COMPOSITE=1.
    mass_comp         = NaN
    t_comp            = 0.0
    comp_status       = "SKIPPED"
    comp_n_iters      = 0
    comp_max_ratio    = NaN
    sl_cmp            = nothing
    sp_cmp            = nothing
    if INCLUDE_COMPOSITE
        t0 = time()
        mass_comp, res_comp, sl_cmp, sp_cmp, comp_n_iters, comp_max_ratio =
            _size_composite(bay, sp_nlp.minimizers)
        t_comp = time() - t0
        comp_status = res_comp === nothing ? "INFEASIBLE" :
                      (res_comp.result_ok ? "OK" : res_comp.solver_status)
        println("  Composite mass  = $(round(mass_comp, digits=1)) kg  ($(round(t_comp, digits=1))s)  " *
                "[$comp_status, iters=$comp_n_iters, max δ/Δ=$(round(comp_max_ratio, digits=2))]")
    end

    mass_ratio = (INCLUDE_COMPOSITE && mass_comp > 0) ? mass_nlp / mass_comp : NaN

    strength_ok  = eval_res.strength_ok
    n_L360_fail  = eval_res.n_L360_fail
    n_L240_fail  = eval_res.n_L240_fail
    n_total      = max(length(sp_nlp.minimizers), 1)
    max_util_M   = eval_res.max_util_M
    max_util_V   = eval_res.max_util_V
    service_ok   = eval_res.serviceability_ok
    max_δ_total  = eval_res.max_δ_total

    println("  NLP-under-composite verdict:")
    println("    strength_ok      = $strength_ok  (max util M/V = $(round(max_util_M, digits=3)) / $(round(max_util_V, digits=3)))")
    println("    serviceability_ok = $service_ok  ($n_L360_fail / $n_total L/360 fail, $n_L240_fail / $n_total L/240 fail)")
    println("    mass_ratio (nlp / comp) = $(round(mass_ratio, digits=3))")

    # Per-beam scatter data (empty when composite was skipped / infeasible).
    beam_rows = (INCLUDE_COMPOSITE && sp_cmp !== nothing) ?
                _per_beam_rows(bay, sl_nlp, sp_nlp, sl_cmp, sp_cmp) :
                NamedTuple[]

    bay_row = (;
        name                  = bay.name,
        scope                 = String(SCOPE),
        slab_type             = String(bay.slab_type),
        slab_sizer            = String(bay.slab_sizer),
        collinear             = bay.collinear,
        max_depth             = bay.max_depth,
        mass_nlp              = mass_nlp,
        mass_comp             = mass_comp,
        comp_status           = comp_status,
        comp_reconcile_iters  = comp_n_iters,
        comp_reconcile_ratio  = comp_max_ratio,
        partial_composite_f   = PARTIAL_COMPOSITE_FACTOR,
        mass_legacy           = ismissing(legacy_mass) ? NaN : legacy_mass,
        mass_ratio            = mass_ratio,
        nlp_strength_ok       = strength_ok,
        nlp_service_ok        = service_ok,
        nlp_max_util_M        = max_util_M,
        nlp_max_util_V        = max_util_V,
        nlp_n_L360_fail       = n_L360_fail,
        nlp_n_L240_fail       = n_L240_fail,
        nlp_max_delta_in      = max_δ_total,
        n_beams               = n_total,
        t_nlp_sec             = t_nlp,
        t_comp_sec            = t_comp,
    )

    return bay_row, beam_rows
end

# ── main ─────────────────────────────────────────────────────────────────────
bays = _scope_bays(SCOPE)
println("Running $(length(bays)) bay(s) under scope = :$SCOPE")
flush(stdout)

# Incremental CSV snapshot so long `:all` sweeps stay monitorable and crash-
# resilient: after each bay we overwrite `per_bay_<scope>.csv` with the rows
# collected so far, and we `flush(stdout)` so `tail -f` on the log file shows
# progress in real time (Julia line-buffers only when stdout is a TTY).
rows           = NamedTuple[]
beam_rows      = NamedTuple[]
csv_out        = joinpath(RESULTS_DIR, "per_bay_$(SCOPE).csv")
csv_beam_out   = joinpath(RESULTS_DIR, "per_beam_$(SCOPE).csv")

for (bi, bay) in enumerate(bays)
    t_bay = time()
    try
        bay_row, beam_row_batch = _compare_one(bay)
        push!(rows, bay_row)
        append!(beam_rows, beam_row_batch)
    catch e
        @warn "Bay failed" bay=bay.name exception=(e, catch_backtrace())
    end
    GC.gc()

    if !isempty(rows)
        CSV.write(csv_out, DataFrame(rows))
    end
    if !isempty(beam_rows)
        CSV.write(csv_beam_out, DataFrame(beam_rows))
    end

    elapsed_min = round((time() - t_bay) / 60, digits=2)
    println("  [$bi / $(length(bays))] $(bay.name) done in $(elapsed_min) min  " *
            "(bays=$(length(rows)), beams=$(length(beam_rows)))")
    flush(stdout)
end

per_bay  = DataFrame(rows)
per_beam = DataFrame(beam_rows)
CSV.write(csv_out,      per_bay)
CSV.write(csv_beam_out, per_beam)
println("\nWrote $csv_out       ($(nrow(per_bay)) rows)")
println("Wrote $csv_beam_out  ($(nrow(per_beam)) rows)")
flush(stdout)

# ── summary table ────────────────────────────────────────────────────────────
if nrow(per_bay) >= 1
    println("\n══════════════════════════════════════════════════════════════")
    println("  Summary across $(nrow(per_bay)) bay(s)")
    println("══════════════════════════════════════════════════════════════")
    finite_ratio = filter(isfinite, per_bay.mass_ratio)
    if !isempty(finite_ratio)
        println("  mass_ratio (NLP / composite):")
        println("    min    = $(round(minimum(finite_ratio), digits=3))")
        println("    median = $(round(median(finite_ratio),  digits=3))")
        println("    mean   = $(round(mean(finite_ratio),    digits=3))")
        println("    max    = $(round(maximum(finite_ratio), digits=3))")
        n_conservative = count(≥(1.0), finite_ratio)
        println("    # bays with mass_nlp ≥ mass_comp: $n_conservative / $(length(finite_ratio))")
    end

    nm = count(per_bay.nlp_strength_ok)
    nd = count(per_bay.nlp_service_ok)
    N  = nrow(per_bay)
    println("  NLP sections under composite FE:")
    println("    strength_ok      = $nm / $N bays")
    println("    serviceability_ok = $nd / $N bays")
end
