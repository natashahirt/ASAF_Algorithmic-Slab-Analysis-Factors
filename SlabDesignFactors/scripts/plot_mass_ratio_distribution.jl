"""
Plot the distribution of NLP/composite mass ratios across all bays in
`per_bay_all.csv` (produced incrementally by
`compare_nlp_vs_composite.jl` under `COMPARE_SCOPE=all`).

Two-panel layout:
  (A) Histogram of `mass_ratio = mass_nlp / mass_comp` with mean,
      median and unity reference.
  (B) `mass_ratio` vs. number of beams per bay — tests whether larger
      bays show a larger composite credit, which is the geometric
      prediction: longer spans → larger Ix demand → larger
      transformed-section contribution.

The script is safe to run while the sweep is still in flight; it just
reports on whatever rows have been written so far.

Usage:
    julia --project=. SlabDesignFactors/scripts/plot_mass_ratio_distribution.jl

Output:
    SlabDesignFactors/plot/figures/composite/mass_ratio_distribution.pdf
"""

include("_scripts.jl")

using CSV, DataFrames, Statistics

CairoMakie.activate!()

const CSV_PATH = joinpath(@__DIR__, "..", "results",
                          "compare_nlp_vs_composite", "per_bay_all.csv")
const FIG_DIR  = joinpath(@__DIR__, "..", "plot", "figures", "composite")
const FIG_PATH = joinpath(FIG_DIR, "mass_ratio_distribution.pdf")

mkpath(FIG_DIR)

isfile(CSV_PATH) || error("Missing $CSV_PATH — has the :all sweep written any rows yet?")

df = CSV.read(CSV_PATH, DataFrame)
df = dropmissing(df, :mass_ratio)
filter!(r -> isfinite(r.mass_ratio), df)

n = nrow(df)
println("Loaded $n bays from $CSV_PATH")
println("  mass_ratio summary:")
println("    min    = $(round(minimum(df.mass_ratio), digits=3))")
println("    p25    = $(round(quantile(df.mass_ratio, 0.25), digits=3))")
println("    median = $(round(median(df.mass_ratio),    digits=3))")
println("    mean   = $(round(mean(df.mass_ratio),      digits=3))")
println("    p75    = $(round(quantile(df.mass_ratio, 0.75), digits=3))")
println("    max    = $(round(maximum(df.mass_ratio), digits=3))")

# ── figure ───────────────────────────────────────────────────────────────────
fontsize       = 14
fontsize_small = 11

fig = Figure(size = (720, 340), fontsize = fontsize)

# ── (A) histogram of mass ratio ──────────────────────────────────────────────
axA = Axis(fig[1, 1],
    xlabel = "NLP / composite steel mass ratio",
    ylabel = "Number of bays",
    title  = "(A) Distribution of paper-vs-composite mass (n = $n bays)",
    xgridvisible = false, ygridvisible = false,
    titlealign = :left,
)

hist!(axA, df.mass_ratio;
      bins        = max(8, min(24, n)),
      color       = (色[:ceruleanblue], 0.75),
      strokewidth = 0.5, strokecolor = :white)

vlines!(axA, [1.0];                 color = :black,     linestyle = :dash,
        linewidth = 1.5, label = "parity (ratio = 1.0)")
vlines!(axA, [median(df.mass_ratio)]; color = :black,   linestyle = :dot,
        linewidth = 1.5,
        label = "median = $(round(median(df.mass_ratio), digits=2))")
vlines!(axA, [mean(df.mass_ratio)];   color = 色[:lilac], linestyle = :solid,
        linewidth = 1.5,
        label = "mean = $(round(mean(df.mass_ratio), digits=2))")

axislegend(axA; position = :rt, framevisible = true, backgroundcolor = :white,
           framecolor = :white, labelsize = fontsize_small,
           patchsize = (14, 10))

# ── (B) mass ratio vs n_beams (bay size proxy) ───────────────────────────────
axB = Axis(fig[1, 2],
    xlabel = "Number of beams per bay",
    ylabel = "NLP / composite mass ratio",
    title  = "(B) Composite credit vs. bay size",
    xgridvisible = false, ygridvisible = false,
    titlealign = :left,
)

scatter!(axB, df.n_beams, df.mass_ratio;
         color      = 色[:ceruleanblue],
         markersize = 8, alpha = 0.8)
hlines!(axB, [1.0]; color = :black, linestyle = :dash, linewidth = 1)

# Annotate the BaU bay if present so readers can locate the r6c4 anchor
bau_idx = findfirst(==("r6c4"), df.name)
if !isnothing(bau_idx)
    bau_row = df[bau_idx, :]
    scatter!(axB, [bau_row.n_beams], [bau_row.mass_ratio];
             color = 色[:lilac], markersize = 12, marker = :diamond,
             label = "BaU (r6c4)")
    text!(axB, bau_row.n_beams, bau_row.mass_ratio;
          text = "  r6c4", align = (:left, :center), fontsize = fontsize_small)
end

save(FIG_PATH, fig)
println("\nWrote figure → $FIG_PATH")
