"""
Replot the strip-resolution convergence figure
(`strip_resolution_convergence.pdf`) at a larger, project-standard font size.

This reads the cached CSV produced by `sensitivity_strip_resolution.jl` and
regenerates the three-panel convergence plot (EC, max deflection, runtime vs
strip spacing) without rerunning the sizing pipeline. Output is written next to
the original under a `_largefont.pdf` suffix so existing figures are preserved.

Fonts are the project defaults (11 / 8) scaled up by ~1.3×:
  * `fontsize        = 14`
  * `fontsize_small  = 11`

Usage (from project root):
    julia --project=. SlabDesignFactors/scripts/replot_strip_resolution.jl
"""

include("_scripts.jl")

CairoMakie.activate!()

# ── paths ─────────────────────────────────────────────────────────────────────
save_path = "SlabDesignFactors/results/sensitivity_strip_resolution/"
csv_path  = joinpath(save_path, "strip_resolution_sensitivity.csv")
out_file  = joinpath(save_path, "strip_resolution_convergence_largefont.pdf")

isfile(csv_path) || error("CSV not found: $csv_path")
results_df = CSV.read(csv_path, DataFrame)

# ── font sizes ────────────────────────────────────────────────────────────────
# Project defaults are 11 / 8; scaled ~1.7× for a print-legible variant where
# the three-panel figure is reproduced at ~half page width.
fontsize       = 18
fontsize_small = 14

# ── runtime outlier flag ──────────────────────────────────────────────────────
# The first two sizing calls of a fresh Julia session pay the full JIT compile
# cost of the Nonconvex/Zygote/Ipopt/JuMP stack (tens of seconds vs <1 s for
# every subsequent call). These show up as the opening rows of the sweep
# (`r1c1 @ spacing=1.0`, `r1c2 @ spacing=1.0`) and otherwise dominate the
# runtime panel. Flag any point whose runtime exceeds 5× the median at the
# same strip spacing so it can be rendered in grey (distinct from the
# colormap-coded converged points). EC and deflection at these rows are
# physically valid, so only the runtime panel is affected.
results_df.runtime_outlier = falses(nrow(results_df))
for s in unique(results_df.spacing)
    idx = findall(==(s), results_df.spacing)
    ts  = results_df.elapsed_s[idx]
    ts_valid = filter(!isnan, ts)
    isempty(ts_valid) && continue
    thresh = 5 * median(ts_valid)
    results_df.runtime_outlier[idx] .= .!isnan.(ts) .& (ts .> thresh)
end
let n_flag = count(results_df.runtime_outlier)
    if n_flag > 0
        flagged = results_df[results_df.runtime_outlier,
                             [:geometry, :spacing, :elapsed_s]]
        @info "Flagging $n_flag runtime outlier(s) as JIT warm-up (grey)" flagged
    end
end

# ── convergence plot (colored by # cells; includes runtime subplot) ───────────
cell_cmap  = Reverse(:dense)   # dark blue = many cells, light blue = few cells
df_plot    = filter(row -> !isnan(row.ec_total), results_df)
cell_range = isempty(df_plot) ? (0, 1) :
             (minimum(df_plot.n_cells), maximum(df_plot.n_cells))

# Match the canvas size of the original convergence figure (and the other
# three-panel figures in SlabDesignFactors/plot/figures/).
fig = Figure(size=(1300, 500), fontsize=fontsize)

axis_kwargs = (
    titlesize      = fontsize,
    xlabelsize     = fontsize,
    ylabelsize     = fontsize,
    xticklabelsize = fontsize_small,
    yticklabelsize = fontsize_small,
)

ax1 = Axis(fig[1, 1];
    xlabel = "Strip spacing (m)",
    ylabel = "Total EC (kgCO₂e/m²)",
    title  = "Embodied carbon vs strip resolution",
    xscale = log10,
    axis_kwargs...,
)
ax2 = Axis(fig[1, 2];
    xlabel = "Strip spacing (m)",
    ylabel = "Max deflection (in)",
    title  = "Max deflection vs strip resolution",
    xscale = log10,
    axis_kwargs...,
)
ax3 = Axis(fig[1, 3];
    xlabel = "Strip spacing (m)",
    ylabel = "Runtime (s)",
    title  = "Runtime vs strip resolution",
    xscale = log10,
    yscale = log10,
    axis_kwargs...,
)

outlier_color = (:grey60, 0.9)

for name in unique(results_df.geometry)
    df_g = filter(row -> row.geometry == name && !isnan(row.ec_total), results_df)
    isempty(df_g) && continue
    sort!(df_g, :spacing)
    n_c = df_g.n_cells[1]

    scatterlines!(ax1, df_g.spacing, df_g.ec_total,
        color=n_c, colormap=cell_cmap, colorrange=cell_range, markersize=8)
    scatterlines!(ax2, df_g.spacing, df_g.max_deflection_in,
        color=n_c, colormap=cell_cmap, colorrange=cell_range, markersize=8)

    # Runtime panel: draw a dotted grey line through all points first, then
    # overlay the solid colour-coded curve through only the converged points.
    # The dotted segment is only visible where it bridges a warm-up point —
    # elsewhere it is painted over by the solid line. This keeps the
    # per-geometry trend readable without hiding the warm-up samples.
    df_ok  = filter(row -> !row.runtime_outlier, df_g)
    df_bad = filter(row ->  row.runtime_outlier, df_g)
    if nrow(df_bad) > 0
        lines!(ax3, df_g.spacing, df_g.elapsed_s;
            color=outlier_color, linestyle=:dot, linewidth=1.5)
    end
    if nrow(df_ok) >= 2
        scatterlines!(ax3, df_ok.spacing, df_ok.elapsed_s,
            color=n_c, colormap=cell_cmap, colorrange=cell_range, markersize=8)
    elseif nrow(df_ok) == 1
        scatter!(ax3, df_ok.spacing, df_ok.elapsed_s,
            color=n_c, colormap=cell_cmap, colorrange=cell_range, markersize=8)
    end
    if nrow(df_bad) > 0
        scatter!(ax3, df_bad.spacing, df_bad.elapsed_s;
            color=outlier_color, markersize=10, marker=:xcross, strokewidth=0)
    end
end

# Grey-marker legend entry for the runtime panel so readers know what the
# off-colour points mean without hunting through the caption.
if any(results_df.runtime_outlier)
    warmup_marker = MarkerElement(color=outlier_color, marker=:xcross, markersize=10)
    axislegend(ax3, [warmup_marker], ["JIT warm-up"];
        position=:rt, framevisible=true, backgroundcolor=:white,
        framecolor=:white, labelsize=fontsize_small, padding=(4, 4, 4, 4))
end

Colorbar(fig[1, 4];
    colormap      = cell_cmap,
    colorrange    = cell_range,
    label         = "# cells in layout",
    labelsize     = fontsize,
    ticklabelsize = fontsize_small,
)

save(out_file, fig)
println("\n✓ Plot saved to $out_file")
display(fig)
