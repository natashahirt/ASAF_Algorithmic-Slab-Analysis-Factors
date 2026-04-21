"""
Per-beam diagnostic for the BaU layout: compare the Nov 2024 bare-steel
NLP (paper approach) against the modern AISC I3.1a transformed-section
composite sizer, and validate the implicit composite-stiffness credit
that the paper's bare-steel L/360 constraint assumes.

For each of the 174 beams in the r6c4 BaU we record, under the actual
tributary load pattern:

  * `A_nlp`, `I_bare_nlp`           — section area and bare steel Ix of
                                       the Nov 2024 minimizer.
  * `I_comp_nlp`                     — transformed-section Ix of the same
                                       minimizer (AISC I3.1a, n=E_s/E_c,
                                       per-side L/8 cap).
  * `A_comp`, `I_bare_comp`, `I_comp_comp` — same three quantities for
                                             the composite-sized
                                             minimizer.

The ratio `I_comp / I_bare` is the *empirical* composite stiffness
credit — i.e. the multiplier the paper's bare-steel L/360 design is
implicitly leaving on the table. If it is ≥ the nominal deflection
reduction factor (DRF) used in the 2024 sweeps, then bare-steel L/360
is a conservative approximation of composite L/360 and the published
embodied-carbon figures are upper bounds.

Outputs
-------
  * `SlabDesignFactors/results/compare_nlp_vs_composite/per_beam_bau.csv`
  * `SlabDesignFactors/plot/figures/composite/drf_validation_bau.pdf`

The figure has two panels:
  (A) per-beam section-area scatter (NLP vs composite, y = x reference)
  (B) histogram of I_comp / I_bare for the NLP minimizers, with the
      DRF = 2.5 reference line marked

Usage (from project root, ~12 min wall time):
    julia --project=. SlabDesignFactors/scripts/paper_drf_validation.jl

Set `PAPER_DRF_FORCE=1` to overwrite an existing CSV cache.
"""

include("_scripts.jl")

using CSV, DataFrames, JSON, Statistics

CairoMakie.activate!()

# ── constants ─────────────────────────────────────────────────────────────────
const GEOM_DIR    = joinpath(@__DIR__, "..", "..", "Geometries", "topology")
const RESULTS_DIR = joinpath(@__DIR__, "..", "results", "compare_nlp_vs_composite")
const FIGURE_DIR  = joinpath(@__DIR__, "..", "plot",  "figures", "composite")
const PER_BEAM_CSV = joinpath(RESULTS_DIR, "per_beam_bau.csv")
const FIGURE_PDF   = joinpath(FIGURE_DIR,  "drf_validation_bau.pdf")

const MAX_DEPTH_IN = 25.0
const BAU = (name="r6c4", slab_type=:isotropic, slab_sizer=:uniform,
             collinear=false, vector_1d=[0.0, 0.0], max_depth=MAX_DEPTH_IN)

mkpath(RESULTS_DIR); mkpath(FIGURE_DIR)

const FORCE_RESIZE = get(ENV, "PAPER_DRF_FORCE", "0") == "1"

# ── helpers shared with compare_nlp_vs_composite.jl ──────────────────────────
function _load_geom(name::AbstractString)
    path = joinpath(GEOM_DIR, name * ".json")
    raw  = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""),
                                  dicttype=Dict))
    return raw isa Dict ? raw : Dict(pairs(raw))
end

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

function _make_sizing_params(; composite::Bool)
    return SlabSizingParams(
        live_load              = psf_to_ksi(50),
        superimposed_dead_load = psf_to_ksi(15),
        slab_dead_load         = 0.0,
        live_factor            = 1.6,
        dead_factor            = 1.2,
        beam_sizer             = :continuous,
        nlp_solver             = :MMA,
        max_depth              = MAX_DEPTH_IN,
        beam_units             = :in,
        serviceability_lim     = 360,
        collinear              = BAU.collinear,
        minimum_continuous     = true,
        n_max_sections         = 0,
        composite_action       = composite,
        partial_composite_factor = 1.0,
        reconcile_max_iter     = 0,
    )
end

"""
    _trib_strips(sizing_params, beam_element) -> (positions, widths, is_perim)

Collect the load-strip positions (0–1) and per-side tributary widths (in)
attached to a beam. Mirrors the logic used inside `postprocess_slab` and
`get_deflection_constraint` so the transformed-section Ix we compute
here matches what the sizer used.
"""
function _trib_strips(sizing_params, beam_element, perim_set)
    loadid_idx = Dict(lid => row for (row, lid) in
                      enumerate(sizing_params.load_df.loadID))
    positions, widths = Float64[], Float64[]
    eid = get_element_id(beam_element)
    for ld in sizing_params.load_dictionary[eid]
        hasproperty(ld, :loadID) || continue
        row = get(loadid_idx, getproperty(ld, :loadID), nothing)
        isnothing(row) && continue
        push!(positions, ld.position)
        push!(widths,    sizing_params.load_df[row, :trib_width])
    end
    idx = findfirst(el -> el === beam_element,
                    sizing_params.model.elements[:beam])
    is_perim = !isnothing(idx) && idx in perim_set
    return positions, widths, is_perim
end

"""
    _per_beam_record(slab_params, sizing_params, label) -> DataFrame

Harvest per-beam section geometry, bare-steel and composite Ix for every
beam in a sized model. Returns a DataFrame with one row per beam.
"""
function _per_beam_record(slab_params, sizing_params, label::AbstractString)
    beam_elements = sizing_params.model.elements[:beam]
    # The bare-steel NLP path leaves `sizing_params.i_perimeter` empty — the
    # perimeter set is only populated inside `optimal_beamsizer` when
    # `composite_action == true`. Fall back to the authoritative
    # `slab_params.i_perimeter` so composite Ix is computed correctly for
    # either regime.
    perim_raw = isempty(sizing_params.i_perimeter) ?
                slab_params.i_perimeter : sizing_params.i_perimeter
    perim_set = Set(perim_raw)
    t_slab_in = maximum(slab_params.slab_depths) *
                convert_to_m[slab_params.slab_units] / convert_to_m[:in]
    E_s = steel_ksi.E
    E_c = sizing_params.E_c

    n = length(beam_elements)
    h_v  = Vector{Float64}(undef, n); w_v  = similar(h_v)
    tw_v = similar(h_v);              tf_v = similar(h_v)
    A_v  = similar(h_v);              Ib_v = similar(h_v);  Ic_v = similar(h_v)
    ratio_v       = similar(h_v)
    n_strips_v    = Vector{Int}(undef, n)
    sum_trib_v    = similar(h_v)
    is_perim_v    = Vector{Bool}(undef, n)
    L_v           = similar(h_v)

    for (i, be) in enumerate(beam_elements)
        vars = sizing_params.minimizers[i]
        h, w, tw, tf = vars
        sec = I_symm(vars...)

        positions, widths, is_perim = _trib_strips(sizing_params, be, perim_set)

        Ix_bare = sec.Ix
        Ix_comp = isempty(widths) ? Ix_bare :
                  get_I_composite_effective(h, w, tw, tf,
                      t_slab_in, E_s, E_c, be.length,
                      positions, widths; is_perimeter=is_perim)

        h_v[i]  = h;  w_v[i]  = w;  tw_v[i] = tw; tf_v[i] = tf
        A_v[i]  = sec.A
        Ib_v[i] = Ix_bare; Ic_v[i] = Ix_comp
        ratio_v[i]    = Ix_comp / Ix_bare
        n_strips_v[i] = length(widths)
        sum_trib_v[i] = sum(widths)
        is_perim_v[i] = is_perim
        L_v[i]        = be.length
    end

    return DataFrame(
        regime            = fill(label, n),
        beam_idx          = collect(1:n),
        L_in              = L_v,
        is_perimeter      = is_perim_v,
        h_in              = h_v,  w_in  = w_v,
        tw_in             = tw_v, tf_in = tf_v,
        A_in2             = A_v,
        Ix_bare_in4       = Ib_v,
        Ix_comp_in4       = Ic_v,
        ratio_Icomp_Ibare = ratio_v,
        n_strips          = n_strips_v,
        sum_trib_in       = sum_trib_v,
    )
end

"""
    _size_bay(; composite, initial_vars) -> (slab_params, sizing_params)

Run either the Nov 2024 bare-steel NLP (`composite=false`) or the modern
transformed-section composite sizer (`composite=true`) on BaU. In
composite mode `initial_vars` should be the NLP minimizers so MMA
warm-starts from the same feasible point that `compare_nlp_vs_composite`
uses; a cold start can park on a much lighter local optimum that fails
the global-FE L/360 check.
"""
function _size_bay(; composite::Bool, initial_vars=Vector{Vector{Float64}}())
    geom_dict = _load_geom(BAU.name)
    slab_params   = _make_slab_params(geom_dict, BAU)
    sizing_params = _make_sizing_params(composite=composite)
    slab_params   = analyze_slab(slab_params)

    if composite
        slab_params, sizing_params =
            optimal_beamsizer(slab_params, sizing_params;
                              initial_vars=initial_vars)
    else
        slab_params, sizing_params = size_bare_nlp!(slab_params, sizing_params)
    end
    return slab_params, sizing_params
end

# ── main: assemble per-beam dataframe ────────────────────────────────────────
if isfile(PER_BEAM_CSV) && !FORCE_RESIZE
    println("Loading cached per-beam data from $PER_BEAM_CSV")
    println("Set PAPER_DRF_FORCE=1 to rerun sizing.")
    combined = CSV.read(PER_BEAM_CSV, DataFrame)
else
    println("[1/2] Sizing BaU with Nov 2024 bare-steel NLP …"); flush(stdout)
    t0 = time()
    sp_nlp, pp_nlp = _size_bay(composite=false)
    println("      done in $(round(time() - t0, digits=1)) s"); flush(stdout)
    df_nlp = _per_beam_record(sp_nlp, pp_nlp, "NLP (bare L/360)")
    println("      assembled per-beam dataframe ($(nrow(df_nlp)) rows)"); flush(stdout)

    println("[2/2] Sizing BaU with AISC I3.1a composite transformed section " *
            "(warm-start from NLP) …"); flush(stdout)
    t0 = time()
    sp_cmp, pp_cmp = _size_bay(composite=true,
                               initial_vars=copy(pp_nlp.minimizers))
    println("      done in $(round(time() - t0, digits=1)) s"); flush(stdout)
    df_cmp = _per_beam_record(sp_cmp, pp_cmp, "Composite (I_eff L/360)")

    combined = vcat(df_nlp, df_cmp)
    CSV.write(PER_BEAM_CSV, combined)
    println("Wrote per-beam data → $PER_BEAM_CSV"); flush(stdout)
end

# ── aggregate statistics reported at stdout ──────────────────────────────────
df_nlp = filter(r -> startswith(r.regime, "NLP"),       combined)
df_cmp = filter(r -> startswith(r.regime, "Composite"), combined)

println("""
\n══════════════════════════════════════════════════════════════
  BaU per-beam comparison  (n = $(nrow(df_nlp)) beams)
══════════════════════════════════════════════════════════════
  Total steel area  [in²]:
    NLP       = $(round(sum(df_nlp.A_in2 .* df_nlp.L_in) * convert_to_m[:in]^3 * ρ_STEEL,   digits=1)) kg
    Composite = $(round(sum(df_cmp.A_in2 .* df_cmp.L_in) * convert_to_m[:in]^3 * ρ_STEEL, digits=1)) kg
  Per-beam I_comp / I_bare (on NLP sections):
    min    = $(round(minimum(df_nlp.ratio_Icomp_Ibare), digits=2))
    p25    = $(round(quantile(df_nlp.ratio_Icomp_Ibare, 0.25), digits=2))
    median = $(round(median(df_nlp.ratio_Icomp_Ibare),  digits=2))
    mean   = $(round(mean(df_nlp.ratio_Icomp_Ibare),    digits=2))
    p75    = $(round(quantile(df_nlp.ratio_Icomp_Ibare, 0.75), digits=2))
    max    = $(round(maximum(df_nlp.ratio_Icomp_Ibare), digits=2))
    fraction ≥ 2.5 : $(round(count(≥(2.5), df_nlp.ratio_Icomp_Ibare) / nrow(df_nlp), digits=3))
""")

# ── figure: 2-panel (area scatter + Icomp/Ibare histogram) ───────────────────
fontsize       = 14
fontsize_small = 11

fig = Figure(size = (720, 340), fontsize=fontsize)

# ── (A) per-beam area scatter NLP vs Composite ───────────────────────────────
axA = Axis(fig[1, 1],
    xlabel = "NLP section area, Aₙₗₚ  [in²]",
    ylabel = "Composite section area, A_cmp  [in²]",
    title  = "(A) Per-beam section area (n = $(nrow(df_nlp)))",
    xgridvisible = false, ygridvisible = false,
    titlealign = :left,
)

a_nlp = df_nlp.A_in2
a_cmp = df_cmp.A_in2
a_max = max(maximum(a_nlp), maximum(a_cmp)) * 1.05
lines!(axA, [0, a_max], [0, a_max],
       color = :black, linestyle = :dash, linewidth = 1)
scatter!(axA, a_nlp, a_cmp;
         color = 色[:ceruleanblue], markersize = 6, alpha = 0.7, transparency = true)
xlims!(axA, 0, a_max); ylims!(axA, 0, a_max)

# annotate mass ratio
m_nlp = sum(df_nlp.A_in2 .* df_nlp.L_in) * convert_to_m[:in]^3 * ρ_STEEL
m_cmp = sum(df_cmp.A_in2 .* df_cmp.L_in) * convert_to_m[:in]^3 * ρ_STEEL
text!(axA, 0.04 * a_max, 0.92 * a_max;
      text  = "mₙₗₚ = $(round(m_nlp, digits=0)) kg\n" *
              "m_cmp = $(round(m_cmp, digits=0)) kg\n" *
              "ratio = $(round(m_nlp / m_cmp, digits=2))",
      align = (:left, :top), fontsize = fontsize_small)

# ── (B) I_comp / I_bare histogram on NLP minimizers ──────────────────────────
axB = Axis(fig[1, 2],
    xlabel = "Empirical stiffness credit, I_comp / I_bare",
    ylabel = "Number of beams",
    title  = "(B) Transformed-section credit on NLP sections",
    xgridvisible = false, ygridvisible = false,
    titlealign = :left,
)

ratios = df_nlp.ratio_Icomp_Ibare
hist!(axB, ratios;
      bins  = 24,
      color = (色[:ceruleanblue], 0.7),
      strokewidth = 0.5, strokecolor = :white)

# DRF = 1.0 (pure bare-steel reference used by the paper's L/360 constraint)
vlines!(axB, [1.0]; color = :black,     linestyle = :solid, linewidth = 1.5,
        label = "DRF = 1.0 (paper)")
# DRF = 2.5 (Nov 2024 sensitivity variant — see `replot_composite_mass_scatter.jl`)
vlines!(axB, [2.5]; color = 色[:lilac], linestyle = :dash, linewidth = 2,
        label = "DRF = 2.5 (sensitivity)")
# Observed central tendency
vlines!(axB, [median(ratios)]; color = :black, linestyle = :dot, linewidth = 1.5,
        label = "median = $(round(median(ratios), digits=2))")

axislegend(axB; position = :rt, framevisible = true, backgroundcolor = :white,
           framecolor = :white, labelsize = fontsize_small,
           patchsize = (14, 10))

save(FIGURE_PDF, fig)
println("\nWrote figure → $FIGURE_PDF")
