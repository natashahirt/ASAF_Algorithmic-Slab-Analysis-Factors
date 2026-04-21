"""
Replot the Monte Carlo best-vs-BaU reduction violins from the CSVs written by
`sensitivity_monte_carlo_ec.jl` (via `dump_reduction_csvs`), without rerunning
the 10 000-sample Monte Carlo loop.

Design
------
For each steel scenario the violin body shows the full reduction distribution.
A central vertical spine runs top-to-bottom across the violin, and a full-
width horizontal median line crosses it left-to-right. Short horizontal dashes
inside the violin mark the interquartile boundaries (P25 and P75). All
numeric annotations (median, P25, P75) are rendered at the same font size,
left-aligned at a fixed x just outside the right edge of the violin, so the
coloured body stays visually uncluttered.

The x-tick for each violin carries the scenario label and the best-layout
identity (e.g. "CLF Typical (P50)\\nbest = grid[6,4]").

A *combined* variant (`plot_reduction_violin_combined` /
`replot_monte_carlo_violin_combined`) overlays multiple design cases
(e.g. `nodeflection_yesslabmin` vs `nodeflection_noslabmin`) as paired
violins grouped by steel scenario, so the headline effect of relaxing the
minimum-slab constraint is visible at a glance.

Usage
-----
Standalone (iterating on the plot after the main MC run has written CSVs):

    julia --project=. SlabDesignFactors/scripts/replot_monte_carlo_violin.jl

Programmatic (called from `sensitivity_monte_carlo_ec.jl`):

    include("replot_monte_carlo_violin.jl")
    replot_monte_carlo_violin(save_path; outpath = ...)
    replot_monte_carlo_violin_combined(save_path_root, design_cases;
                                       outpath = ...)
"""

# Pull in the shared plotting stack / palette (`色`) only when invoked
# standalone; when `include`d from the MC script these are already loaded.
if !@isdefined(CairoMakie)
    include("_scripts.jl")
end
using Random: MersenneTwister
CairoMakie.activate!()

const _VIOLIN_FONTSIZE      = 14
const _VIOLIN_SMALLFONTSIZE = 11
# Inside-violin P25 / median / P75 annotations — set smaller than the axis /
# tick labels so the coloured bodies stay visually uncluttered.
const _VIOLIN_ANNOT_FONTSIZE = 11
const _VIOLIN_PALETTE = [色[:ceruleanblue], 色[:charcoalgrey], 色[:magenta]]
# Two-tone palette for the combined design-case plot; chosen to read as
# "baseline vs relaxed" rather than tracking scenario identity.
const _VIOLIN_DESIGN_PALETTE = [色[:ceruleanblue], 色[:magenta]]

"""
    plot_reduction_violin_from_data(samples_df, summary_df; outpath,
                                    palette = _VIOLIN_PALETTE,
                                    title_suffix = "")

Render the annotated violin plot from precomputed samples + summary data.

`samples_df` is wide (one column per scenario name, one row per MC draw).
`summary_df` has one row per scenario, in display order, with columns
`scenario`, `label`, `best_rowcol`, `p5`, `p25`, `p50`, `p75`, `p95`.
`palette` is indexed in the same order as `summary_df` rows.
"""
function plot_reduction_violin_from_data(samples_df::DataFrame,
                                         summary_df::DataFrame;
                                         outpath::AbstractString,
                                         palette = _VIOLIN_PALETTE,
                                         title_suffix::AbstractString = "")
    n = nrow(summary_df)

    tick_labels = ["$(row.label)\nbest = $(row.best_rowcol)"
                   for row in eachrow(summary_df)]

    fig = Figure(size = (820, 500), fontsize = _VIOLIN_FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "Steel scenario",
        ylabel = "EC reduction best-vs-BaU (%)",
        title  = "Reduction distribution across steel scenarios$title_suffix",
        xticks = (1:n, tick_labels),
        xticklabelsize     = _VIOLIN_SMALLFONTSIZE,
        xlabelpadding      = 18,
        topspinevisible    = false,
        rightspinevisible  = false,
        bottomspinevisible = false,
        xticksvisible      = false,
    )

    violin_width = 1.0
    quartile_hw  = 0.05   # em-dash length for P25 / P75 notches

    all_vals = Float64[]

    for (k, row) in enumerate(eachrow(summary_df))
        col     = palette[mod1(k, length(palette))]
        samples = samples_df[!, row.scenario]
        append!(all_vals, samples)

        vmin, vmax = extrema(samples)

        # `show_median = true` lets CairoMakie fit the median line to the
        # violin body width at p50, so it never spills outside the KDE.
        violin!(ax, fill(k, length(samples)), samples,
            color       = (col, 0.45),
            strokecolor = col, strokewidth = 1.2,
            width       = violin_width,
            mediancolor = col,
            show_median = true)

        # Central spine — top-to-bottom across the full violin body.
        lines!(ax, [k, k], [vmin, vmax],
            color = col, linewidth = 2)

        # P25 / P75 notches on the spine.
        for q in (row.p25, row.p75)
            lines!(ax, [k - quartile_hw, k + quartile_hw], [q, q],
                color = col, linewidth = 2)
        end
    end

    pad = 0.08 * (maximum(all_vals) - minimum(all_vals))
    ylims!(ax, (minimum(all_vals) - pad, maximum(all_vals) + pad))
    xlims!(ax, (0.4, n + 0.6))

    save(outpath, fig)
    return fig
end

"""
    replot_monte_carlo_violin(save_path; outpath = ..., kwargs...)

Convenience entry point: reads `reduction_samples.csv` and
`reduction_summary.csv` from `save_path` and renders the violin to `outpath`
(default: `save_path/violin_reduction_by_steel_scenario.pdf`).
"""
function replot_monte_carlo_violin(save_path::AbstractString;
                                   outpath::AbstractString = joinpath(save_path,
                                       "violin_reduction_by_steel_scenario.pdf"),
                                   kwargs...)
    samples_df = CSV.read(joinpath(save_path, "reduction_samples.csv"), DataFrame)
    summary_df = CSV.read(joinpath(save_path, "reduction_summary.csv"), DataFrame)
    title_suffix = "  (N=$(nrow(samples_df)))"
    return plot_reduction_violin_from_data(samples_df, summary_df;
        outpath = outpath, title_suffix = title_suffix, kwargs...)
end

"""
    load_layout_reductions(save_path, summary_df) -> Dict{String, Vector{Float64}}

For each scenario in `summary_df`, compute the *cross-layout* reduction vs BaU
from the on-disk per-scenario results CSV plus `scenario_summary.csv`.

This is the quantity the layout-space plot visualises: one reduction value per
valid design in the search space, capturing the design-choice variance that
the ECC-only MC violin cannot express.
"""
function load_layout_reductions(save_path::AbstractString, summary_df::DataFrame)
    scen_summary = CSV.read(joinpath(save_path, "scenario_summary.csv"), DataFrame)
    bau_by_scen  = Dict(row.scenario => row.bau_ec_nominal for row in eachrow(scen_summary))

    out = Dict{String, Vector{Float64}}()
    for row in eachrow(summary_df)
        results_path = joinpath(save_path, row.scenario, "monte_carlo_ec_results.csv")
        results = CSV.read(results_path, DataFrame)
        bau_ec  = bau_by_scen[row.scenario]
        out[row.scenario] = (bau_ec .- results.ec_nominal) ./ bau_ec .* 100
    end
    return out
end

"""
    plot_reduction_layoutspace_from_data(layout_reductions, samples_df, summary_df;
                                         outpath, palette, title_suffix,
                                         show_points, point_alpha)

Render a *layout-space* reduction violin per steel scenario, with the MC 90%
CI of the best-vs-BaU reduction overlaid as a narrow rangebar on the right
flank. The figure surfaces two distinct sources of uncertainty that the plain
MC violin conflates / hides:

  • Design-choice variance — fat violin body + jittered layout dots, built
    from one nominal reduction-vs-BaU per valid layout in the search space.
  • ECC-coefficient variance on the nominal best — narrow rangebar + gold
    star at the nominal reduction, drawn from the same MC samples the
    original violin used.

`layout_reductions` is a Dict keyed by scenario name to a Vector of per-layout
reductions-vs-BaU (%). `samples_df` and `summary_df` share the schemas
consumed by `plot_reduction_violin_from_data`.
"""
function plot_reduction_layoutspace_from_data(layout_reductions::AbstractDict,
                                              samples_df::DataFrame,
                                              summary_df::DataFrame;
                                              outpath::AbstractString,
                                              palette = _VIOLIN_PALETTE,
                                              title_suffix::AbstractString = "",
                                              show_points::Bool = true,
                                              point_alpha::Real = 0.35)
    # `samples_df` is accepted (and ignored) so the call site mirrors
    # `plot_reduction_violin_from_data`; the layout-space view no longer
    # overlays the ECC MC envelope since it's an order of magnitude tighter
    # than the design-space distribution and only clutters the figure.
    _ = samples_df
    n = nrow(summary_df)

    tick_labels = ["$(row.label)\nbest = $(row.best_rowcol)"
                   for row in eachrow(summary_df)]

    fig = Figure(size = (900, 520), fontsize = _VIOLIN_FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "Steel scenario",
        ylabel = "EC reduction vs BaU (%)",
        title  = "Layout-space reduction distribution$title_suffix",
        xticks = (1:n, tick_labels),
        xticklabelsize     = _VIOLIN_SMALLFONTSIZE,
        xlabelpadding      = 18,
        topspinevisible    = false,
        rightspinevisible  = false,
        bottomspinevisible = false,
        xticksvisible      = false,
    )

    # Zero-reduction reference so "beats BaU" is visually obvious.
    hlines!(ax, [0.0], color = :black, linewidth = 1)

    violin_width = 0.80
    jitter_w     = 0.18
    quartile_hw  = 0.05   # em-dash-length notches for P25 / P75

    all_vals = Float64[]

    for (k, row) in enumerate(eachrow(summary_df))
        col   = palette[mod1(k, length(palette))]
        lvals = layout_reductions[row.scenario]
        append!(all_vals, lvals)

        vmin, vmax = extrema(lvals)

        # Violin body. `show_median = true` lets CairoMakie fit the median
        # line to the body width at p50 so it never escapes the violin.
        violin!(ax, fill(k, length(lvals)), lvals,
            color       = (col, 0.22),
            strokecolor = col, strokewidth = 1.0,
            width       = violin_width,
            mediancolor = col,
            show_median = true)

        # Jittered layout dots — readers can literally count designs beating BaU.
        if show_points
            rng = MersenneTwister(4242 + k)
            xs  = k .+ (rand(rng, length(lvals)) .- 0.5) .* jitter_w
            scatter!(ax, xs, lvals,
                color       = (col, point_alpha),
                strokewidth = 0,
                markersize  = 3)
        end

        # Central spine — top-to-bottom across the full violin body.
        lines!(ax, [k, k], [vmin, vmax],
            color = col, linewidth = 2)

        # P25 / P75 notches on the spine.
        for q in quantile(lvals, (0.25, 0.75))
            lines!(ax, [k - quartile_hw, k + quartile_hw], [q, q],
                color = col, linewidth = 2)
        end
    end

    pad = 0.08 * (maximum(all_vals) - minimum(all_vals))
    ylims!(ax, (minimum(all_vals) - pad, maximum(all_vals) + pad))
    xlims!(ax, (0.4, n + 0.6))

    save(outpath, fig)
    return fig
end

"""
    replot_monte_carlo_layoutspace(save_path; outpath = ..., kwargs...)

Convenience entry point for the layout-space plot. Reads
`reduction_samples.csv`, `reduction_summary.csv`, `scenario_summary.csv`, and
the per-scenario `<scenario>/monte_carlo_ec_results.csv` under `save_path`;
assembles the per-layout reductions; and renders the figure to `outpath`
(default: `save_path/layoutspace_reduction_by_steel_scenario.pdf`).
"""
function replot_monte_carlo_layoutspace(save_path::AbstractString;
                                        outpath::AbstractString = joinpath(save_path,
                                            "layoutspace_reduction_by_steel_scenario.pdf"),
                                        kwargs...)
    samples_df = CSV.read(joinpath(save_path, "reduction_samples.csv"), DataFrame)
    summary_df = CSV.read(joinpath(save_path, "reduction_summary.csv"), DataFrame)
    layout_red = load_layout_reductions(save_path, summary_df)
    n_layouts  = length(first(values(layout_red)))
    title_suffix = "  (N_MC=$(nrow(samples_df)), N_layouts=$n_layouts)"
    return plot_reduction_layoutspace_from_data(layout_red, samples_df, summary_df;
        outpath = outpath, title_suffix = title_suffix, kwargs...)
end

"""
    plot_reduction_violin_combined(design_cases_data; outpath,
                                   palette = _VIOLIN_DESIGN_PALETTE,
                                   title_suffix = "")

Paired-violin variant comparing two (or more) design cases side-by-side within
each steel scenario group. Mirrors the single-case vocabulary (violin body,
central spine, P25 / P75 notches, in-body median stroke); numeric labels and
best-layout tags are omitted to keep paired bodies legible. The design-case
legend sits underneath the axis as a horizontal strip.

`design_cases_data` is a Vector of NamedTuples
`(name, label, samples_df, summary_df)`, in display order (e.g. with-slab-min
first, no-slab-min second). All entries must share the same ordered set of
scenarios (by `summary_df.scenario`); an error is raised otherwise.
"""
function plot_reduction_violin_combined(design_cases_data::AbstractVector;
                                        outpath::AbstractString,
                                        palette = _VIOLIN_DESIGN_PALETTE,
                                        title_suffix::AbstractString = "")
    @assert !isempty(design_cases_data) "need at least one design case"
    n_designs = length(design_cases_data)
    scenarios = design_cases_data[1].summary_df.scenario
    for dc in design_cases_data
        @assert dc.summary_df.scenario == scenarios "all design cases must share scenarios in matching order"
    end
    n_scen      = length(scenarios)
    scen_labels = design_cases_data[1].summary_df.label

    # Position pairs within each scenario group; paired violins sit on either
    # side of the group center with a small interior gap.
    pair_half = 0.24
    offsets = n_designs == 1 ? [0.0] :
              n_designs == 2 ? [-pair_half, pair_half] :
              collect(range(-pair_half, pair_half, length = n_designs))

    # Narrower bodies than the single-case plot so two fit per group.
    violin_width = 0.38
    quartile_hw  = 0.04

    fig = Figure(size = (1050, 620), fontsize = _VIOLIN_FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "Steel scenario",
        ylabel = "EC reduction best-vs-BaU (%)",
        title  = "Reduction distribution: with vs without slab-min constraint$title_suffix",
        xticks = (1:n_scen, scen_labels),
        xticklabelsize     = _VIOLIN_SMALLFONTSIZE,
        xlabelpadding      = 18,
        topspinevisible    = false,
        rightspinevisible  = false,
        bottomspinevisible = false,
        xticksvisible      = false,
    )

    all_vals = Float64[]

    for (d, dc) in enumerate(design_cases_data)
        col = palette[mod1(d, length(palette))]

        for (g, _) in enumerate(scenarios)
            samples = dc.samples_df[!, dc.summary_df.scenario[g]]
            append!(all_vals, samples)
            x          = g + offsets[d]
            vmin, vmax = extrema(samples)
            srow       = dc.summary_df[g, :]

            violin!(ax, fill(x, length(samples)), samples,
                color       = (col, 0.45),
                strokecolor = col, strokewidth = 1.2,
                width       = violin_width,
                mediancolor = col,
                show_median = true)

            lines!(ax, [x, x], [vmin, vmax], color = col, linewidth = 2)
            for q in (srow.p25, srow.p75)
                lines!(ax, [x - quartile_hw, x + quartile_hw], [q, q],
                    color = col, linewidth = 2)
            end
        end
    end

    pad = 0.12 * (maximum(all_vals) - minimum(all_vals))
    ylims!(ax, (minimum(all_vals) - pad, maximum(all_vals) + pad))
    xlims!(ax, (0.3, n_scen + 0.7))

    # Horizontal legend underneath the axis.
    legend_elements = [PolyElement(
            color       = (palette[mod1(d, length(palette))], 0.45),
            strokecolor =  palette[mod1(d, length(palette))],
            strokewidth = 1.2,
        ) for d in 1:n_designs]
    legend_labels = [dc.label for dc in design_cases_data]
    Legend(fig[2, 1], legend_elements, legend_labels;
        orientation    = :horizontal,
        framevisible   = false,
        labelsize      = _VIOLIN_SMALLFONTSIZE,
        tellheight     = true,
        tellwidth      = false,
        patchsize      = (18, 12),
        colgap         = 28,
        padding        = (0, 0, 0, 8),
    )

    save(outpath, fig)
    return fig
end

"""
    replot_monte_carlo_violin_combined(save_path_root, design_cases;
                                       outpath = ..., kwargs...)

Convenience entry point: reads `reduction_samples.csv` and
`reduction_summary.csv` from each `save_path_root/<design_case.name>/` and
renders the combined violin to `outpath`.

`design_cases` is a Vector of NamedTuples with `name` and `label` fields,
matching the list used in `sensitivity_monte_carlo_ec.jl`.
"""
function replot_monte_carlo_violin_combined(save_path_root::AbstractString,
                                            design_cases;
                                            outpath::AbstractString = joinpath(save_path_root,
                                                "violin_reduction_combined.pdf"),
                                            kwargs...)
    data = map(design_cases) do dc
        dir        = joinpath(save_path_root, dc.name)
        samples_df = CSV.read(joinpath(dir, "reduction_samples.csv"), DataFrame)
        summary_df = CSV.read(joinpath(dir, "reduction_summary.csv"), DataFrame)
        (name = dc.name, label = dc.label, samples_df = samples_df, summary_df = summary_df)
    end
    n_samples    = nrow(data[1].samples_df)
    title_suffix = "  (N=$n_samples per violin)"
    return plot_reduction_violin_combined(data;
        outpath = outpath, title_suffix = title_suffix, kwargs...)
end

"""
    plot_reduction_layoutspace_combined(design_cases_data; outpath, palette,
                                        title_suffix)

Paired layout-space violins, one pair per steel scenario, with the MC 90% CI
of each design case's best-vs-BaU rendered as a small rangebar on the inner
flank of its violin. Intended for the `nodeflection_yesslabmin` vs
`nodeflection_noslabmin` comparison so relaxing the slab-min constraint is
visible both as a shift in *what designs are feasible* (violin shape) and as
a shift in *what the best-layout achieves* (rangebar + star).

`design_cases_data` is a Vector of NamedTuples
`(name, label, layout_reductions, samples_df, summary_df)` in display order.
All entries must share the same ordered set of scenarios.
"""
function plot_reduction_layoutspace_combined(design_cases_data::AbstractVector;
                                             outpath::AbstractString,
                                             palette = _VIOLIN_DESIGN_PALETTE,
                                             title_suffix::AbstractString = "",
                                             show_points::Bool = true,
                                             point_alpha::Real = 0.30)
    @assert !isempty(design_cases_data) "need at least one design case"
    n_designs = length(design_cases_data)
    scenarios = design_cases_data[1].summary_df.scenario
    for dc in design_cases_data
        @assert dc.summary_df.scenario == scenarios "all design cases must share scenarios in matching order"
    end
    n_scen      = length(scenarios)
    scen_labels = design_cases_data[1].summary_df.label

    pair_half = 0.24
    offsets = n_designs == 1 ? [0.0] :
              n_designs == 2 ? [-pair_half, pair_half] :
              collect(range(-pair_half, pair_half, length = n_designs))

    violin_width = 0.38
    jitter_w     = 0.08
    quartile_hw  = 0.04   # em-dash-length notches on the central spine

    _ = title_suffix  # title retired; suffix kept in API for back-compat.

    fig = Figure(size = (1050, 600), fontsize = _VIOLIN_FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "Steel scenario",
        ylabel = "EC reduction vs BaU (%)",
        xticks = (1:n_scen, scen_labels),
        xticklabelsize     = _VIOLIN_SMALLFONTSIZE,
        xlabelpadding      = 18,
        topspinevisible    = false,
        rightspinevisible  = false,
        bottomspinevisible = false,
        xticksvisible      = false,
    )
    hlines!(ax, [0.0], color = :black, linewidth = 1)

    all_vals = Float64[]

    for (d, dc) in enumerate(design_cases_data)
        col = palette[mod1(d, length(palette))]

        for (g, scen) in enumerate(scenarios)
            lvals = dc.layout_reductions[scen]
            append!(all_vals, lvals)

            x = g + offsets[d]
            vmin, vmax = extrema(lvals)

            # Violin body; `show_median = true` guarantees the median line
            # tracks the body width at p50 (it cannot spill outside the KDE).
            violin!(ax, fill(x, length(lvals)), lvals,
                color       = (col, 0.22),
                strokecolor = col, strokewidth = 1.0,
                width       = violin_width,
                mediancolor = col,
                show_median = true)

            # Jittered layout dots — narrower jitter than the single-case
            # version so paired violins don't overlap.
            if show_points
                rng = MersenneTwister(4242 + 100*d + g)
                xs  = x .+ (rand(rng, length(lvals)) .- 0.5) .* jitter_w
                scatter!(ax, xs, lvals,
                    color       = (col, point_alpha),
                    strokewidth = 0,
                    markersize  = 2.5)
            end

            # Central spine — top-to-bottom across the full violin body.
            lines!(ax, [x, x], [vmin, vmax],
                color = col, linewidth = 1.6)

            # P25 / P75 notches on the spine.
            for q in quantile(lvals, (0.25, 0.75))
                lines!(ax, [x - quartile_hw, x + quartile_hw], [q, q],
                    color = col, linewidth = 1.6)
            end
        end
    end

    pad = 0.08 * (maximum(all_vals) - minimum(all_vals))
    ylims!(ax, (minimum(all_vals) - pad, maximum(all_vals) + pad))
    xlims!(ax, (0.3, n_scen + 0.7))

    # Only the design-case legend remains — the two colors need decoding.
    design_elements = [PolyElement(
            color       = (palette[mod1(d, length(palette))], 0.22),
            strokecolor =  palette[mod1(d, length(palette))],
            strokewidth = 1.0,
        ) for d in 1:n_designs]
    design_labels = [dc.label for dc in design_cases_data]
    Legend(fig[2, 1], design_elements, design_labels;
        orientation  = :horizontal,
        framevisible = false,
        labelsize    = _VIOLIN_SMALLFONTSIZE,
        tellheight   = true,
        tellwidth    = false,
        patchsize    = (18, 12),
        colgap       = 28,
        padding      = (0, 0, 0, 8),
    )

    save(outpath, fig)
    return fig
end

"""
    replot_monte_carlo_layoutspace_combined(save_path_root, design_cases;
                                            outpath = ..., kwargs...)

Convenience entry point: loads per-design-case CSVs (samples, summary,
per-scenario results) from `save_path_root/<design_case.name>/` and renders
the combined layout-space plot to `outpath`.
"""
function replot_monte_carlo_layoutspace_combined(save_path_root::AbstractString,
                                                 design_cases;
                                                 outpath::AbstractString = joinpath(save_path_root,
                                                     "layoutspace_reduction_combined.pdf"),
                                                 kwargs...)
    data = map(design_cases) do dc
        dir        = joinpath(save_path_root, dc.name)
        samples_df = CSV.read(joinpath(dir, "reduction_samples.csv"), DataFrame)
        summary_df = CSV.read(joinpath(dir, "reduction_summary.csv"), DataFrame)
        layout_red = load_layout_reductions(dir, summary_df)
        (name              = dc.name,
         label             = dc.label,
         layout_reductions = layout_red,
         samples_df        = samples_df,
         summary_df        = summary_df)
    end
    n_samples    = nrow(data[1].samples_df)
    n_layouts    = length(first(values(data[1].layout_reductions)))
    title_suffix = "  (N_MC=$n_samples, N_layouts=$n_layouts per violin)"
    return plot_reduction_layoutspace_combined(data;
        outpath = outpath, title_suffix = title_suffix, kwargs...)
end

# Standalone entry point — only fires when the file is executed directly,
# not when it is `include`d from `sensitivity_monte_carlo_ec.jl`. Regenerates
# every per-design-case violin plus the combined comparison from the CSVs
# already on disk.
if abspath(PROGRAM_FILE) == @__FILE__
    save_path_root = "SlabDesignFactors/results/sensitivity_monte_carlo_ec/"

    # Kept in sync with `design_cases` in `sensitivity_monte_carlo_ec.jl`.
    design_cases = [
        (name = "nodeflection_yesslabmin", label = "Slab minimum of 0.125m"),
        (name = "nodeflection_noslabmin",  label = "No slab minimum"),
    ]

    for dc in design_cases
        dir        = joinpath(save_path_root, dc.name)
        violin_out = joinpath(dir, "violin_reduction_by_steel_scenario.pdf")
        layout_out = joinpath(dir, "layoutspace_reduction_by_steel_scenario.pdf")
        replot_monte_carlo_violin(dir; outpath = violin_out)
        println("✓ Per-case violin regenerated → $violin_out")
        replot_monte_carlo_layoutspace(dir; outpath = layout_out)
        println("✓ Per-case layout-space plot regenerated → $layout_out")
    end

    combined_violin = joinpath(save_path_root, "violin_reduction_combined.pdf")
    replot_monte_carlo_violin_combined(save_path_root, design_cases;
        outpath = combined_violin)
    println("✓ Combined design-case violin regenerated → $combined_violin")

    combined_layout = joinpath(save_path_root, "layoutspace_reduction_combined.pdf")
    replot_monte_carlo_layoutspace_combined(save_path_root, design_cases;
        outpath = combined_layout)
    println("✓ Combined design-case layout-space plot regenerated → $combined_layout")
end
