"""
Total-mass scatter + paper-BaU benchmark for the topology slab set.

For each of the four large-scale design sweeps (the paper's "remote results"
CSVs), compare the per-topology steel mass against three current-pipeline
variants (naive-2.5×, full-2.5×, composite). Reports regression stats
(slope-through-origin, R², median y/x, MAPE) so the baseline whose sizing
is closest to the composite reference can be identified — this is the
baseline the paper should adopt for its BaU definition.

Series definitions
------------------
* **Legacy baselines** — `steel_norm × area` from the four
  `remote_results_<def><smin>/topology.csv` CSVs, each filtered to the
  isotropic + uniform + discrete + collinear + depth = 40 in. configuration:
      yesdef_yessmin : `remote_results_yesdeflection_yesslabmin`
      yesdef_nosmin  : `remote_results_yesdeflection_noslabmin`
      nodef_yessmin  : `remote_results_nodeflection_yesslabmin`
      nodef_nosmin   : `remote_results_nodeflection_noslabmin`
* **Naive 2.5×** — current pipeline's re-derivation of the pre-composite
  sizing method, i.e. bare steel with the deflection constraint *active* and
  `deflection_reduction_factor = 2.5` baked into the sizer
  (`(δ_FE + δ_SW)/2.5 ≤ L/360` during optimization). This is the sizing
  approach that existed in the codebase before composite action was
  introduced.
* **Full 2.5×** — identical sizing method as *Naive 2.5×*, but the column
  values come from the published paper baseline CSV via
  `replot_composite_mass_scatter.jl`'s `mass_bare25` cache. Kept as a
  reference so drift between the current pipeline and the paper run is
  visible in one figure.
* **Composite** — current pipeline, `composite_action=true` (cache's
  `mass_composite`) — the reference against which each legacy baseline is
  benchmarked.

Caching
-------
Reuses `composite_mass_cache.csv` (depth = 40, isotropic, uniform, discrete,
collinear) for the bare/full25/composite columns. The naive-2.5× column is
sized fresh on first run and persisted to `naive25_mass_cache_v2.csv` (the
`_v2` suffix avoids reusing the strength-only cache from the prior
definition). Set `REPLOT_FORCE_RESIZE = 1` to force a re-run.

Usage (from project root):
    julia --project=. SlabDesignFactors/scripts/compare_4way_mass_scatter.jl
"""

include("_scripts.jl")
using PrettyTables

CairoMakie.activate!()

# ── plotting style (mirrors replot_composite_mass_scatter.jl) ────────────────
fontsize       = 14
fontsize_small = 11

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

const SCATTER_KWARGS = (
    markersize   = 4,
    alpha        = 0.6,
    transparency = true,
)

# ── paths ────────────────────────────────────────────────────────────────────
main_path         = "Geometries/topology/"
save_path         = "SlabDesignFactors/plot/figures/composite/"
cache_path        = "SlabDesignFactors/results/composite_mass_cache.csv"
# Cache file is versioned (`_v2`) because the `naive25` variant now sizes with
# the deflection constraint active (pre-composite method). The v1 cache held
# strength-only masses and would otherwise be silently reused.
naive_cache_path  = "SlabDesignFactors/results/naive25_mass_cache_v2.csv"

# Legacy paper-baseline CSVs. Each one is the output of a full design-sweep
# run under a different (deflection × slab-min) regime; we load them all so
# the regression-vs-composite benchmark can tell us which regime is closest
# to the composite reference. Keys are the short names used in column labels,
# filenames, and the stats table.
const BASELINES = [
    ("yesdef_yessmin",
     "SlabDesignFactors/results/remote_results_yesdeflection_yesslabmin/topology.csv",
     "Yes-defl • Yes-slab-min"),
    ("yesdef_nosmin",
     "SlabDesignFactors/results/remote_results_yesdeflection_noslabmin/topology.csv",
     "Yes-defl • No-slab-min"),
    ("nodef_yessmin",
     "SlabDesignFactors/results/remote_results_nodeflection_yesslabmin/topology.csv",
     "No-defl • Yes-slab-min"),
    ("nodef_nosmin",
     "SlabDesignFactors/results/remote_results_nodeflection_noslabmin/topology.csv",
     "No-defl • No-slab-min"),
]

mkpath(save_path)
out_file       = joinpath(save_path, "composite_mass_scatter_4way.pdf")
out_benchmark  = joinpath(save_path, "baseline_vs_composite_benchmark.pdf")
stats_out_path = "SlabDesignFactors/results/comparison_composite/baseline_vs_composite_stats.csv"

force_resize = get(ENV, "REPLOT_FORCE_RESIZE", "0") == "1"

# ── load existing cache (mass_bare, mass_bare25, mass_composite) ─────────────
isfile(cache_path) || error(
    "Composite cache not found at $cache_path — run replot_composite_mass_scatter.jl first.",
)
summary_df = CSV.read(cache_path, DataFrame)
println("Loaded composite cache: $(nrow(summary_df)) geometries.")

# ── load each legacy baseline CSV → `mass_paper_<key>` column ────────────────
#
# Every CSV is filtered to the same canonical topology configuration
# (isotropic / uniform / discrete / collinear / vector=(0,0) / depth=40 in)
# so the rows align geometry-by-geometry with the composite cache.
"""
    load_baseline!(summary_df, csv_path, key)

Add `mass_paper_<key>` to `summary_df` by joining `csv_path` on the
`geometry` column, filtering to the canonical configuration, and computing
`steel_norm × area` (kg) per geometry. Leaves NaN for missing rows.
"""
function load_baseline!(summary_df::DataFrame, csv_path::AbstractString,
                        key::AbstractString)
    isfile(csv_path) || error("Baseline CSV not found at $csv_path.")
    df = CSV.read(csv_path, DataFrame)
    subset = filter(row ->
        row.slab_type   == "isotropic" &&
        row.slab_sizer  == "uniform"   &&
        row.beam_sizer  == "discrete"  &&
        row.collinear   == true        &&
        row.vector_1d_x == 0.0         &&
        row.vector_1d_y == 0.0         &&
        row.max_depth   == 40.0        &&
        row.area > 0,
        df,
    )
    col = Symbol("mass_paper_" * key)
    summary_df[!, col] = fill(NaN, nrow(summary_df))
    for i in 1:nrow(summary_df)
        m = filter(r -> r.name == summary_df.geometry[i], subset)
        if nrow(m) == 1 && m[1, :steel_norm] > 0
            summary_df[i, col] = m[1, :steel_norm] * m[1, :area]
        end
    end
    return count(!isnan, summary_df[!, col])
end

for (key, path, label) in BASELINES
    n = load_baseline!(summary_df, path, key)
    println("Loaded baseline '$key' [$label]: $n / $(nrow(summary_df)) geometries.")
end

# The first baseline (yesdef_yessmin) is the "paper" series used by the
# legacy 3-panel scatter figure below — alias for backwards compatibility.
summary_df.mass_paper = summary_df[!, Symbol("mass_paper_" * BASELINES[1][1])]

# ── compute / load naive-2.5× column (pre-composite sizing method) ───────────
# Reproduces the sizer approach used before composite action was introduced:
# bare steel (`composite_action=false`) with the deflection constraint active
# and `drf=2.5` baked into the sizer, so the optimization enforces
# `(δ_FE + δ_SW)/2.5 ≤ L/360` rather than the strict bare-section L/360.
# Geometry/load/material settings mirror the composite cache so only the
# sizing method differs between the scatter panels.
function size_naive25_for_geometry(name::AbstractString, geometry_dict)
    geom_c, _ = Base.invokelatest(generate_from_json, geometry_dict;
                                  plot=false, drawn=false)
    slab_params = SlabAnalysisParams(
        geom_c,
        slab_name     = String(name),
        slab_type     = :isotropic,
        vector_1d     = [1.0, 0.0],  # ignored for :isotropic; kept aligned with the composite cache
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
        composite_action            = false,
        deflection_reduction_factor = 2.5,
        deflection_limit            = true,
    )
    slab_params = analyze_slab(slab_params)
    slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)
    isempty(sizing_params.minimizers) && return 0.0
    scaled_beam_lengths = [be.length for be in sizing_params.model.elements[:beam]]
    return sum(I_symm(m...).A * scaled_beam_lengths[i]
               for (i, m) in enumerate(sizing_params.minimizers)
              ) * convert_to_m[:in]^3 * ρ_STEEL
end

naive_use_cache = isfile(naive_cache_path) && !force_resize
if naive_use_cache
    println("Using cached naive-2.5× sizing from $naive_cache_path")
    println("Set REPLOT_FORCE_RESIZE=1 to force a fresh sizing run.")
    naive_df = CSV.read(naive_cache_path, DataFrame)
    summary_df.mass_naive25 = fill(NaN, nrow(summary_df))
    for i in 1:nrow(summary_df)
        m = filter(r -> r.geometry == summary_df.geometry[i], naive_df)
        if nrow(m) == 1
            summary_df.mass_naive25[i] = m[1, :mass_naive25]
        end
    end
else
    println("\nSizing naive-2.5× (pre-composite method, δ/2.5 ≤ L/360 active) for $(nrow(summary_df)) geometries…")
    summary_df.mass_naive25 = zeros(Float64, nrow(summary_df))
    for (i, name) in enumerate(summary_df.geometry)
        path = joinpath(main_path, "$(name).json")
        if !isfile(path)
            @warn "Missing geometry JSON; leaving mass_naive25 = 0" name
            continue
        end
        raw = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""),
                                    dicttype=Dict))
        geometry_dict = raw isa Dict ? raw : Dict(pairs(raw))
        t0 = time()
        try
            summary_df.mass_naive25[i] = size_naive25_for_geometry(name, geometry_dict)
        catch e
            @warn "Naive-2.5× sizing failed; leaving mass=0" name exception=e
        end
        dt = time() - t0
        println("  [$i/$(nrow(summary_df))] $name  mass=" *
                "$(round(summary_df.mass_naive25[i], digits=1)) kg  ($(round(dt, digits=1))s)")
        GC.gc()
    end
    CSV.write(naive_cache_path,
              DataFrame(geometry = summary_df.geometry,
                        mass_naive25 = summary_df.mass_naive25))
    println("Cached naive-2.5× results to $naive_cache_path")
end

# ── build 3-panel scatter: paper vs each new-pipeline variant ────────────────
if nrow(summary_df) >= 2
    L_tot = summary_df.total_beam_length
    len_max = max(maximum(L_tot), 1e-9)
    len_range = (0.0, len_max)
    len_cmap = :blues

    panels = [
        (:mass_naive25,   "Naive 2.5× (pre-composite sizer)"),
        (:mass_bare25,    "Full 2.5× (sized with δ-limit)"),
        (:mass_composite, "Composite"),
    ]

    fig = Figure(size=(190 * 6, 190 * 2.4), fontsize=fontsize)

    for (pi, (col_y, label_y)) in enumerate(panels)
        mx = summary_df.mass_paper
        my = summary_df[!, col_y]
        # Both axes must be > 0 for a comparison: a zero on either side means
        # the corresponding sizing failed for that geometry.
        nz = findall(i -> !isnan(mx[i]) && mx[i] > 0 && my[i] > 0, eachindex(mx))
        length(nz) < 2 && continue

        mx_nz = mx[nz]
        my_nz = my[nz]
        L_nz  = L_tot[nz]

        ax = Axis(fig[1, pi];
            xlabel         = "Mass — Original (paper) (kg)",
            ylabel         = "Mass — $label_y (kg)",
            titlesize      = fontsize,
            xlabelsize     = fontsize,     ylabelsize     = fontsize,
            xticklabelsize = fontsize_small, yticklabelsize = fontsize_small,
            aspect         = AxisAspect(1),
        )
        scatter!(ax, mx_nz, my_nz;
                 color=L_nz, colormap=len_cmap, colorrange=len_range,
                 SCATTER_KWARGS...)

        # Origin-anchored best fit: my ≈ slope · mx
        slope = sum(mx_nz .* my_nz) / sum(mx_nz .* mx_nz)
        lim = max(maximum(mx_nz), maximum(my_nz)) * 1.05
        xl = range(0.0, lim; length=50)
        fit_line = lines!(ax, xl, slope .* xl;
                          color=色[:magenta], linewidth=1.5)
        ref_line = lines!(ax, [0, lim], [0, lim];
                          color=(:black, 0.4), linestyle=:dash, linewidth=1)
        scatter_marker = MarkerElement(color=色[:ceruleanblue], marker=:circle,
                                       markersize=fontsize_small)
        axislegend(ax,
                   [scatter_marker, fit_line, ref_line],
                   ["Geometries (n=$(length(nz)))",
                    "Best fit (slope=$(round(slope, digits=3)))",
                    "1:1"];
                   LEGEND_KWARGS...)
    end

    Colorbar(fig[1, length(panels) + 1];
        colormap      = len_cmap,
        colorrange    = len_range,
        label         = "Total beam length (m)",
        labelsize     = fontsize,
        ticklabelsize = fontsize_small,
    )
    Label(fig[0, :], "Design mass: paper baseline vs current-pipeline variants";
          fontsize=fontsize, font=:bold, halign=:center)
    save(out_file, fig)
    println("\nSaved: $out_file")
else
    @warn "Not enough geometries to plot (need ≥ 2)."
end

# ── Benchmark: which legacy baseline is closest to composite? ────────────────
#
# For each baseline (x) vs. each current-pipeline method (y), compute
# regression-through-origin statistics so we can populate a Table-12-style
# summary. The statistic we care about for the paper BaU choice is the
# row with method = Composite: whichever baseline's slope is closest to
# 1.0 (with the lowest MAPE) is the regime whose topology.csv most
# faithfully reproduces the composite reference mass.
"""
    benchmark_stats(x, y)

Return a NamedTuple of regression statistics for `y ~ slope · x` fit
through the origin, along with per-point percent-error metrics (MAPE
and median y/x). Skips NaN / non-positive entries.
"""
function benchmark_stats(x::AbstractVector, y::AbstractVector)
    valid = [(x[i], y[i]) for i in eachindex(x)
             if !isnan(x[i]) && !isnan(y[i]) && x[i] > 0 && y[i] > 0]
    n = length(valid)
    n < 2 && return (; n, slope=NaN, r2=NaN, median_yx=NaN, mape=NaN)

    xv = [p[1] for p in valid]
    yv = [p[2] for p in valid]
    slope      = sum(xv .* yv) / sum(xv .^ 2)
    ss_res     = sum((yv .- slope .* xv) .^ 2)
    ss_tot     = sum((yv .- mean(yv)) .^ 2)
    r2         = ss_tot > 0 ? 1 - ss_res / ss_tot : NaN
    ratios     = yv ./ xv
    median_yx  = median(ratios)
    mape       = mean(abs.(ratios .- 1.0)) * 100
    return (; n, slope, r2, median_yx, mape)
end

METHODS = [
    (:mass_naive25,   "Naive 2.5×"),
    (:mass_bare25,    "Full 2.5×"),
    (:mass_composite, "Composite"),
]

stats_rows = DataFrame(
    baseline  = String[],
    label     = String[],
    method    = String[],
    n         = Int[],
    slope     = Float64[],
    r2        = Float64[],
    median_yx = Float64[],
    mape_pct  = Float64[],
)

for (key, _, label) in BASELINES, (method_col, method_label) in METHODS
    x = summary_df[!, Symbol("mass_paper_" * key)]
    y = summary_df[!, method_col]
    s = benchmark_stats(x, y)
    push!(stats_rows, (key, label, method_label, s.n,
                       round(s.slope,     digits=3),
                       round(s.r2,        digits=3),
                       round(s.median_yx, digits=3),
                       round(s.mape,      digits=1)))
end

println("\n══════════════════════════════════════════════")
println("  Paper-baseline regression against current-pipeline methods")
println("  (slope = y/x through origin; 1.0 is a perfect match.")
println("   MAPE computed in the same reference frame.)")
println("══════════════════════════════════════════════\n")
pretty_table(stats_rows;
             column_labels = ["Baseline", "Label", "Method", "n",
                              "Slope", "R²", "Median y/x", "MAPE %"],
             alignment = [:l, :l, :l, :r, :r, :r, :r, :r])

# Headline: for each baseline, the Composite-row tells us how close that
# legacy sweep is to the composite reference. Print a compact summary.
println("\n── Composite-vs-baseline headline (lower MAPE + slope near 1.0 = best) ──")
comp_rows = filter(r -> r.method == "Composite", stats_rows)
sort!(comp_rows, :mape_pct)
for r in eachrow(comp_rows)
    flag = r === first(eachrow(comp_rows)) ? "  ← closest to composite" : ""
    println("  $(rpad(r.label, 26))  slope=$(r.slope)  R²=$(r.r2)  " *
            "MAPE=$(r.mape_pct)%$flag")
end

CSV.write(stats_out_path, stats_rows)
println("\nStats CSV → $stats_out_path")

# ── Compact benchmark figure: Composite vs. each legacy baseline ─────────────
#
# One panel per baseline, composite on the y-axis. Scatter points are
# colored by total beam length as in the 3-panel figure above. The fit
# slope and MAPE are annotated in the legend so the reader can pick the
# baseline at a glance.
if nrow(summary_df) >= 2
    L_tot     = summary_df.total_beam_length
    len_max   = max(maximum(L_tot), 1e-9)
    len_range = (0.0, len_max)
    len_cmap  = :blues

    fig_b = Figure(size=(190 * length(BASELINES) * 1.2, 190 * 2.6),
                   fontsize=fontsize)

    for (pi, (key, _, label)) in enumerate(BASELINES)
        mx = summary_df[!, Symbol("mass_paper_" * key)]
        my = summary_df.mass_composite
        nz = findall(i -> !isnan(mx[i]) && mx[i] > 0 && my[i] > 0,
                     eachindex(mx))
        length(nz) < 2 && continue

        mx_nz, my_nz, L_nz = mx[nz], my[nz], L_tot[nz]
        s = benchmark_stats(mx_nz, my_nz)

        ax = Axis(fig_b[1, pi];
            xlabel         = "Baseline mass (kg) — $label",
            ylabel         = "Composite mass (kg)",
            titlesize      = fontsize,
            xlabelsize     = fontsize,     ylabelsize     = fontsize,
            xticklabelsize = fontsize_small, yticklabelsize = fontsize_small,
            aspect         = AxisAspect(1),
        )
        scatter!(ax, mx_nz, my_nz;
                 color=L_nz, colormap=len_cmap, colorrange=len_range,
                 SCATTER_KWARGS...)

        lim      = max(maximum(mx_nz), maximum(my_nz)) * 1.05
        xl       = range(0.0, lim; length=50)
        fit_line = lines!(ax, xl, s.slope .* xl;
                          color=色[:magenta], linewidth=1.5)
        ref_line = lines!(ax, [0, lim], [0, lim];
                          color=(:black, 0.4), linestyle=:dash, linewidth=1)

        scatter_marker = MarkerElement(color=色[:ceruleanblue], marker=:circle,
                                       markersize=fontsize_small)
        axislegend(ax,
                   [scatter_marker, fit_line, ref_line],
                   ["n = $(length(nz))",
                    "slope = $(round(s.slope, digits=3))   " *
                    "MAPE = $(round(s.mape, digits=1))%",
                    "1:1"];
                   LEGEND_KWARGS...)
    end

    Colorbar(fig_b[1, length(BASELINES) + 1];
        colormap      = len_cmap,
        colorrange    = len_range,
        label         = "Total beam length (m)",
        labelsize     = fontsize,
        ticklabelsize = fontsize_small,
    )
    Label(fig_b[0, :],
          "Composite reference vs. each legacy paper baseline";
          fontsize=fontsize, font=:bold, halign=:center)
    save(out_benchmark, fig_b)
    println("\nSaved: $out_benchmark")
end

println("\nDone.")
