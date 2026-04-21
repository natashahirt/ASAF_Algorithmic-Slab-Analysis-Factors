"""
Per-beam mass scatter across **all** bays of the `:all` sweep, once
`compare_nlp_vs_composite.jl` has emitted `per_beam_all.csv`.

Each point is one beam; x = NLP bare-steel mass, y = composite mass.
Extends the BaU-only `plot_perbeam_mass_bau.jl` figure to the full
topology set (n = 36 bays × ~100-250 beams each).

Two panels:
  (A) All beams on a single scatter with a y = x parity reference,
      colour-coded by topology to catch bay-specific clustering.
  (B) Distribution of per-beam `m_nlp - m_cmp` (kg saved per beam by
      composite sizing), shown as a histogram with zero/median
      references — complements the total-mass view by revealing how
      concentrated the savings are in girder-class beams.

Safe to run on a partial CSV (the compare script writes
`per_beam_all.csv` incrementally after every bay).

Usage:
    julia --project=. SlabDesignFactors/scripts/plot_perbeam_mass_all.jl

Output:
    SlabDesignFactors/plot/figures/composite/perbeam_mass_all.pdf
"""

include("_scripts.jl")

using CSV, DataFrames, Statistics

CairoMakie.activate!()

const CSV_PATH = joinpath(@__DIR__, "..", "results",
                          "compare_nlp_vs_composite", "per_beam_all.csv")
const FIG_DIR  = joinpath(@__DIR__, "..", "plot", "figures", "composite")
const FIG_PATH = joinpath(FIG_DIR, "perbeam_mass_all.pdf")

mkpath(FIG_DIR)

isfile(CSV_PATH) || error("Missing $CSV_PATH — rerun compare_nlp_vs_composite.jl " *
                          "under COMPARE_SCOPE=all with the updated script that " *
                          "emits per-beam data.")

df = CSV.read(CSV_PATH, DataFrame)
df = dropmissing(df, [:m_nlp_kg, :m_cmp_kg])
filter!(r -> isfinite(r.m_nlp_kg) && isfinite(r.m_cmp_kg), df)

n_beams = nrow(df)
bays    = unique(df.name)
n_bays  = length(bays)

println("Loaded $(n_beams) beams across $(n_bays) bays from $CSV_PATH")
println("  Totals   - NLP: $(round(sum(df.m_nlp_kg), digits=0)) kg, " *
        "Composite: $(round(sum(df.m_cmp_kg), digits=0)) kg, " *
        "ratio: $(round(sum(df.m_nlp_kg) / sum(df.m_cmp_kg), digits=3))")
Δ = df.m_nlp_kg .- df.m_cmp_kg
println("  Per-beam savings (m_nlp - m_cmp):")
for (k, v) in (("min", minimum(Δ)), ("p25", quantile(Δ, 0.25)),
               ("median", median(Δ)),
               ("p75", quantile(Δ, 0.75)), ("max", maximum(Δ)))
    println("    $k  = $(round(v, digits=2)) kg")
end
println("  fraction of beams with m_nlp == m_cmp : " *
        "$(round(count(iszero, Δ) / length(Δ), digits=3))")

# ── figure ───────────────────────────────────────────────────────────────────
fontsize       = 14
fontsize_small = 11

fig = Figure(size = (900, 440), fontsize = fontsize)

# ── (A) per-beam scatter, coloured by bay ────────────────────────────────────
axA = Axis(fig[1, 1],
    xlabel = "Nov 2024 NLP beam mass [kg]",
    ylabel = "Composite beam mass [kg]",
    title  = "(A) Per-beam mass, all topologies (n = $(n_beams) beams, $(n_bays) bays)",
    xgridvisible = false, ygridvisible = false,
    titlealign = :left,
    aspect = AxisAspect(1),
)

m_max = max(maximum(df.m_nlp_kg), maximum(df.m_cmp_kg)) * 1.05
lines!(axA, [0, m_max], [0, m_max];
       color = :black, linestyle = :dash, linewidth = 1)

# Distinct hue per bay via a cycled colormap. The viridis family reads
# well under b&w printing and avoids the issues with categorical
# palettes when n_bays is large (~36 here).
bay_ids      = Dict(b => i for (i, b) in enumerate(bays))
colour_grad  = cgrad(:viridis, n_bays; categorical = true)
bay_colour_v = [colour_grad[bay_ids[name]] for name in df.name]

scatter!(axA, df.m_nlp_kg, df.m_cmp_kg;
         color = bay_colour_v, markersize = 5, alpha = 0.55,
         transparency = true, strokewidth = 0)

xlims!(axA, 0, m_max); ylims!(axA, 0, m_max)

text!(axA, 0.04 * m_max, 0.94 * m_max;
      text  = "∑mₙₗₚ = $(round(sum(df.m_nlp_kg) / 1000, digits=1)) t\n" *
              "∑m_cmp = $(round(sum(df.m_cmp_kg) / 1000, digits=1)) t\n" *
              "Σ ratio = $(round(sum(df.m_nlp_kg) / sum(df.m_cmp_kg), digits=2))",
      align = (:left, :top), fontsize = fontsize_small)

# ── (B) histogram of per-beam composite savings ──────────────────────────────
axB = Axis(fig[1, 2],
    xlabel = "Per-beam composite savings, mₙₗₚ − m_cmp [kg]",
    ylabel = "Number of beams",
    title  = "(B) Where the composite credit lands",
    xgridvisible = false, ygridvisible = false,
    titlealign = :left,
)

hist!(axB, Δ; bins = 40, color = (色[:ceruleanblue], 0.75),
      strokewidth = 0.5, strokecolor = :white)
vlines!(axB, [0.0]; color = :black, linestyle = :dash, linewidth = 1.5,
        label = "no change")
vlines!(axB, [median(Δ)]; color = :black, linestyle = :dot, linewidth = 1.5,
        label = "median = $(round(median(Δ), digits=1)) kg")
vlines!(axB, [mean(Δ)];   color = 色[:lilac], linestyle = :solid,
        linewidth = 1.5,
        label = "mean = $(round(mean(Δ), digits=1)) kg")

axislegend(axB; position = :rt, framevisible = true, backgroundcolor = :white,
           framecolor = :white, labelsize = fontsize_small,
           patchsize = (14, 10))

save(FIG_PATH, fig)
println("\nWrote figure → $FIG_PATH")
