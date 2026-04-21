"""
Bay-level mass scatter: total composite steel mass (kg) vs. total Nov
2024 NLP steel mass (kg), one point per topology from the `:all`
sweep's incremental `per_bay_all.csv`.

Replicates the per-beam BaU plot (`plot_perbeam_mass_bau.jl`) at the
bay level. Points above the parity line represent topologies where the
bare-steel paper approach is heavier than the modern composite sizer,
i.e. where the L/360 bare-steel assumption is conservative.

Safe to run while the sweep is still in flight; it simply reports on
whatever rows have been written so far. Set `ANNOTATE_NAMES=1` to
label every bay with its `<row><col>` ID.

Usage:
    julia --project=. SlabDesignFactors/scripts/plot_perbay_mass_scatter.jl

Output:
    SlabDesignFactors/plot/figures/composite/perbay_mass_scatter.pdf
"""

include("_scripts.jl")

using CSV, DataFrames, Statistics

CairoMakie.activate!()

const CSV_PATH = joinpath(@__DIR__, "..", "results",
                          "compare_nlp_vs_composite", "per_bay_all.csv")
const FIG_DIR  = joinpath(@__DIR__, "..", "plot", "figures", "composite")
const FIG_PATH = joinpath(FIG_DIR, "perbay_mass_scatter.pdf")

const ANNOTATE_NAMES = get(ENV, "ANNOTATE_NAMES", "0") == "1"

mkpath(FIG_DIR)

isfile(CSV_PATH) || error("Missing $CSV_PATH — has the :all sweep written any rows yet?")

df = CSV.read(CSV_PATH, DataFrame)
df = dropmissing(df, [:mass_nlp, :mass_comp])
filter!(r -> isfinite(r.mass_nlp) && isfinite(r.mass_comp) && r.mass_comp > 0, df)
sort!(df, :name)

n = nrow(df)
sum_nlp = sum(df.mass_nlp)
sum_cmp = sum(df.mass_comp)

println("Loaded $n bays from $CSV_PATH")
println("  total — NLP: $(round(sum_nlp, digits=0)) kg,  " *
        "Composite: $(round(sum_cmp, digits=0)) kg,  " *
        "ratio: $(round(sum_nlp / sum_cmp, digits=3))")
println("  mass_ratio summary (nlp / comp):")
for (k, v) in (("min", minimum(df.mass_ratio)),
               ("p25", quantile(df.mass_ratio, 0.25)),
               ("median", median(df.mass_ratio)),
               ("mean", mean(df.mass_ratio)),
               ("p75", quantile(df.mass_ratio, 0.75)),
               ("max", maximum(df.mass_ratio)))
    println("    $k  = $(round(v, digits=3))")
end

# ── figure ───────────────────────────────────────────────────────────────────
fontsize       = 14
fontsize_small = 11

fig = Figure(size = (460, 430), fontsize = fontsize)

axA = Axis(fig[1, 1],
    xlabel = "Nov 2024 NLP bay mass [kg]",
    ylabel = "Composite bay mass [kg]",
    title  = "Per-bay mass, all topologies (n = $n)",
    xgridvisible = false, ygridvisible = false,
    titlealign = :left,
    aspect = AxisAspect(1),
)

m_max = max(maximum(df.mass_nlp), maximum(df.mass_comp)) * 1.05
lines!(axA, [0, m_max], [0, m_max];
       color = :black, linestyle = :dash, linewidth = 1, label = "y = x (parity)")

# Highlight the BaU bay if present so readers can locate r6c4 as the anchor.
is_bau = df.name .== "r6c4"
scatter!(axA, df.mass_nlp[.!is_bau], df.mass_comp[.!is_bau];
         color = 色[:ceruleanblue], markersize = 8, alpha = 0.85)
if any(is_bau)
    scatter!(axA, df.mass_nlp[is_bau], df.mass_comp[is_bau];
             color = 色[:lilac], markersize = 12, marker = :diamond)
end

if ANNOTATE_NAMES
    for r in eachrow(df)
        text!(axA, r.mass_nlp, r.mass_comp;
              text = "  " * r.name, align = (:left, :center),
              fontsize = fontsize_small - 2, color = (:black, 0.6))
    end
end

xlims!(axA, 0, m_max); ylims!(axA, 0, m_max)

text!(axA, 0.04 * m_max, 0.94 * m_max;
      text  = "∑mₙₗₚ = $(round(sum_nlp, digits=0)) kg\n" *
              "∑m_cmp = $(round(sum_cmp, digits=0)) kg\n" *
              "Σ ratio = $(round(sum_nlp / sum_cmp, digits=2))\n" *
              "median ratio = $(round(median(df.mass_ratio), digits=2))",
      align = (:left, :top), fontsize = fontsize_small)

elem_bay = MarkerElement(color = 色[:ceruleanblue], marker = :circle,  markersize = 8)
elem_bau = MarkerElement(color = 色[:lilac],       marker = :diamond, markersize = 10)
elem_par = LineElement(color = :black, linestyle = :dash, linewidth = 1)

legend_entries = any(is_bau) ? [elem_bay, elem_bau, elem_par] :
                               [elem_bay, elem_par]
legend_labels  = any(is_bau) ? ["topology bay", "BaU (r6c4)", "y = x (parity)"] :
                               ["topology bay", "y = x (parity)"]
Legend(fig[2, 1], legend_entries, legend_labels;
       orientation = :horizontal, tellheight = true, tellwidth = false,
       framevisible = false, labelsize = fontsize_small, patchsize = (14, 10))

save(FIG_PATH, fig)
println("\nWrote figure → $FIG_PATH")
