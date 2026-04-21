"""
Replot the composite total-mass scatter (`composite_mass_scatter.pdf`) at a
larger, project-standard font size and with the same plotting conventions as
the other SlabDesignFactors figures (Arial theme, project legend/scatter
styling, shared `色` colour dict).

This is a lean variant of `compare_composite.jl` that:
  * reruns the beam-sizing pipeline for all topology geometries under the same
    three cases (bare steel, bare steel 2.5×, composite),
  * skips the deflection-profile computation and per-geometry plots,
  * writes only the total-mass scatter to a new filename so existing figures
    in `SlabDesignFactors/plot/figures/composite/` are preserved.

Fonts are the project defaults (11 / 8) scaled up by ~1.3×:
  * `fontsize        = 14`
  * `fontsize_small  = 11`

Usage (from project root):
    julia --project=. SlabDesignFactors/scripts/replot_composite_mass_scatter.jl
"""

include("_scripts.jl")

CairoMakie.activate!()
# Arial theme is installed globally by SlabDesignFactors/plot/plotting/utils.jl,
# which is loaded via `_scripts.jl` → `_plotting.jl`. No extra set_theme! here.

# ── cases ─────────────────────────────────────────────────────────────────────
# Colours come from the shared `色` dict (defined in plot/plotting/utils.jl).
cases = [
    (label="Bare steel",        color=色[:charcoalgrey], composite=false, drf=1.0),
    (label="Bare steel (2.5×)", color=色[:lilac],        composite=false, drf=2.5),
    (label="Composite",         color=色[:ceruleanblue], composite=true,  drf=1.0),
]

# Project defaults are 11 / 8; scaled up ~1.3× for the larger-font variant.
fontsize       = 14
fontsize_small = 11

# Shared axislegend kwargs that mirror the other SlabDesignFactors plots.
const LEGEND_KWARGS = (
    position        = :lt,
    orientation     = :vertical,
    labelhalign     = :left,
    framevisible    = true,
    backgroundcolor = :white,
    framecolor      = :white,
    labelsize       = fontsize_small,
    patchsize       = (2, 10),
    padding         = (0, 0, 0, 0),
)

# Scatter kwargs matching the project "markersize ~4, alpha, transparency" style.
const SCATTER_KWARGS = (
    markersize   = 4,
    alpha        = 0.6,
    transparency = true,
)

# ── paths ─────────────────────────────────────────────────────────────────────
main_path  = "Geometries/topology/"
save_path  = "SlabDesignFactors/plot/figures/composite/"
cache_path = "SlabDesignFactors/results/composite_mass_cache.csv"
# Original (paper) DRF=2.5 baseline. Used to overwrite `mass_bare25` in the
# cache so the scatter compares against the published bare-2.5× result rather
# than the current pipeline's re-derivation.
old_baseline_path = "SlabDesignFactors/results/remote_results_yesdeflection_yesslabmin/topology.csv"
mkpath(save_path)
mkpath(dirname(cache_path))
out_file   = joinpath(save_path, "composite_mass_scatter_largefont.pdf")

# Skip the sizing pipeline if a cached summary already exists. Set the env var
# `REPLOT_FORCE_RESIZE=1` to force a fresh run.
use_cache = isfile(cache_path) && get(ENV, "REPLOT_FORCE_RESIZE", "0") != "1"

summary_df = DataFrame(
    geometry          = String[],
    total_beam_length = Float64[],
    mass_bare         = Float64[],
    mass_bare25       = Float64[],
    mass_composite    = Float64[],
)

if use_cache
    println("Using cached sizing results from $cache_path")
    println("Set REPLOT_FORCE_RESIZE=1 to force a fresh sizing run.")
    summary_df = CSV.read(cache_path, DataFrame)
else

json_files = filter(x -> endswith(x, ".json"), readdir(main_path))
isempty(json_files) && @warn "No JSON files found in $main_path"

# ── collect per-geometry totals for each case ────────────────────────────────
for json_file in json_files
    name = replace(json_file, ".json" => "")
    path = joinpath(main_path, json_file)

    println("\n════════════════════════════════════════════")
    println("  $name")
    println("════════════════════════════════════════════")

    raw = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))
    geometry_dict = raw isa Dict ? raw : Dict(pairs(raw))

    # Build geometry once to record total beam length; each case rebuilds its
    # own model so sizing is independent.
    geom_ref, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)
    total_L = sum(be.length for be in geom_ref.elements[:beam])

    case_masses = zeros(Float64, length(cases))

    for (ci, c) in enumerate(cases)
        geom_c, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)

        slab_params = SlabAnalysisParams(
            geom_c,
            slab_name     = name,
            slab_type     = :isotropic,
            vector_1d     = [1.0, 0.0],
            slab_sizer    = :uniform,
            spacing       = 0.1,
            plot_analysis = false,
            fix_param     = true,
            slab_units    = :m,
        )

        sizing_params = SlabSizingParams(
            live_load                   = psf_to_ksi(50),
            superimposed_dead_load      = psf_to_ksi(15),
            slab_dead_load              = 0.0,
            live_factor                 = 1.6,
            dead_factor                 = 1.2,
            beam_sizer                  = :discrete,
            max_depth                   = 40.0,
            beam_units                  = :in,
            serviceability_lim          = 360,
            collinear                   = true,
            minimum_continuous          = true,
            n_max_sections              = 0,
            composite_action            = c.composite,
            deflection_reduction_factor = c.drf,
        )

        slab_params = analyze_slab(slab_params)

        t0 = time()
        slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)
        dt = time() - t0

        if !isempty(sizing_params.minimizers)
            scaled_beam_lengths = [be.length for be in sizing_params.model.elements[:beam]]
            case_masses[ci] = sum(
                I_symm(m...).A * scaled_beam_lengths[i]
                for (i, m) in enumerate(sizing_params.minimizers)
            ) * convert_to_m[:in]^3 * ρ_STEEL
        end

        println("  [$(c.label)] sized in $(round(dt, digits=2))s  " *
                "mass=$(round(case_masses[ci], digits=1)) kg")
    end

    push!(summary_df, (name, total_L,
                       case_masses[1], case_masses[2], case_masses[3]))

    GC.gc()
end

CSV.write(cache_path, summary_df)
println("Cached sizing results to $cache_path")

end  # if !use_cache

# ── Override `mass_bare25` with the published paper baseline ─────────────────
# The cache's bare-2.5× column is a re-derivation by the *current* pipeline,
# which has drifted from the paper's run. For figure consistency, substitute
# the equivalent total mass (= steel_norm × area) from the original
# `remote_results_yesdeflection_yesslabmin/topology.csv` baseline, filtered to
# the closest matching configuration: isotropic + uniform + discrete +
# collinear at 40 in. depth (vector_1d=[0,0]).
if isfile(old_baseline_path)
    old_df = CSV.read(old_baseline_path, DataFrame)
    old_subset = filter(row ->
        row.slab_type   == "isotropic" &&
        row.slab_sizer  == "uniform"   &&
        row.beam_sizer  == "discrete"  &&
        row.collinear   == true        &&
        row.vector_1d_x == 0.0         &&
        row.vector_1d_y == 0.0         &&
        row.max_depth   == 40.0        &&
        row.area > 0,
        old_df,
    )
    global n_overridden = 0
    for i in 1:nrow(summary_df)
        m = filter(r -> r.name == summary_df.geometry[i], old_subset)
        if nrow(m) == 1 && m[1, :steel_norm] > 0
            summary_df.mass_bare25[i] = m[1, :steel_norm] * m[1, :area]
            global n_overridden += 1
        end
    end
    println("Overrode mass_bare25 from old baseline for ",
            "$n_overridden / $(nrow(summary_df)) geometries.")
else
    @warn "Old baseline CSV not found; keeping pipeline-derived mass_bare25." path=old_baseline_path
end

# ══════════════════════════════════════════════════════════════════════════════
#  TOTAL-MASS SCATTER (three panels)
# ══════════════════════════════════════════════════════════════════════════════

mass_pairs = [
    (1, 2, "Bare steel",        "Bare steel (2.5×)"),
    (1, 3, "Bare steel",        "Composite"),
    (2, 3, "Bare steel (2.5×)", "Composite"),
]

if nrow(summary_df) >= 2
    L_tot = summary_df.total_beam_length
    len_max = max(maximum(L_tot), 1e-9)
    len_range = (0.0, len_max)
    len_cmap = :blues
    m_cols = [:mass_bare, :mass_bare25, :mass_composite]

    # Panel size follows the ~190 px module used elsewhere in the project,
    # scaled up so the three panels + colorbar read comfortably at this font.
    fig_mass = Figure(size=(190 * 6, 190 * 2.4), fontsize=fontsize)

    for (pi, (ci_x, ci_y, lab_x, lab_y)) in enumerate(mass_pairs)
        mx = summary_df[!, m_cols[ci_x]]
        my = summary_df[!, m_cols[ci_y]]
        # Require BOTH cases to have a successful sizing — a zero on either axis
        # means that case failed for that geometry and would otherwise drag the
        # best-fit line toward the origin.
        nz_m = findall(i -> mx[i] > 0 && my[i] > 0, eachindex(mx))
        length(nz_m) < 2 && continue

        mx_nz = mx[nz_m]
        my_nz = my[nz_m]
        L_nz  = L_tot[nz_m]

        ax = Axis(fig_mass[1, pi];
            xlabel         = "Mass — $lab_x (kg)",
            ylabel         = "Mass — $lab_y (kg)",
            titlesize      = fontsize,
            xlabelsize     = fontsize,     ylabelsize     = fontsize,
            xticklabelsize = fontsize_small, yticklabelsize = fontsize_small,
            aspect         = AxisAspect(1),
        )
        scatter!(ax, mx_nz, my_nz;
                 color=L_nz, colormap=len_cmap, colorrange=len_range,
                 SCATTER_KWARGS...)

        # Best-fit line through the origin (no intercept). For mass-vs-mass
        # comparisons this is the physically meaningful trend: when one design
        # has zero mass, the other should too. The slope is the least-squares
        # estimator of `my ≈ slope · mx`, i.e. slope = (mxᵀ my) / (mxᵀ mx).
        slope = sum(mx_nz .* my_nz) / sum(mx_nz .* mx_nz)
        lim = max(maximum(mx_nz), maximum(my_nz)) * 1.05
        xl = range(0.0, lim; length=50)
        yl = slope .* xl
        fit_line = lines!(ax, xl, yl;
                          color=色[:magenta], linewidth=1.5,
                          label="Best fit (slope=$(round(slope, digits=3)))")
        ref_line = lines!(ax, [0, lim], [0, lim];
                          color=(:black, 0.4), linestyle=:dash, linewidth=1,
                          label="1:1")
        scatter_marker = MarkerElement(color=色[:ceruleanblue], marker=:circle,
                                       markersize=fontsize_small)
        axislegend(ax,
                   [scatter_marker, fit_line, ref_line],
                   ["Geometries", "Best fit (slope=$(round(slope, digits=3)))", "1:1"];
                   LEGEND_KWARGS...)
    end

    Colorbar(fig_mass[1, 4];
        colormap      = len_cmap,
        colorrange    = len_range,
        label         = "Total beam length (m)",
        labelsize     = fontsize,
        ticklabelsize = fontsize_small,
    )
    Label(fig_mass[0, :], "Design mass comparison (total)";
          fontsize=fontsize, font=:bold, halign=:center)
    save(out_file, fig_mass)
    println("\nSaved: $out_file")
else
    @warn "Not enough geometries ran successfully to produce scatter (need ≥ 2)."
end

println("\nDone.")
