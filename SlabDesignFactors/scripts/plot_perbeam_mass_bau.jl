"""
Per-beam mass scatter: composite steel mass (kg) vs. NLP steel mass (kg)
for every beam in the BaU layout. Reads the cached per-beam dataframe
`per_beam_bau.csv` produced by `paper_drf_validation.jl` and emits a
single-panel figure for the paper.

Each point is one beam; the dashed `y = x` line marks parity, and
points below the line are beams where the composite sizer found a
lighter feasible section than the Nov 2024 bare-steel NLP.

Usage:
    julia --project=. SlabDesignFactors/scripts/plot_perbeam_mass_bau.jl

Output:
    SlabDesignFactors/plot/figures/composite/perbeam_mass_bau.pdf
"""

include("_scripts.jl")

using CSV, DataFrames, Statistics

CairoMakie.activate!()

const CSV_PATH = joinpath(@__DIR__, "..", "results",
                          "compare_nlp_vs_composite", "per_beam_bau.csv")
const FIG_DIR  = joinpath(@__DIR__, "..", "plot", "figures", "composite")
const FIG_PATH = joinpath(FIG_DIR, "perbeam_mass_bau.pdf")

mkpath(FIG_DIR)

isfile(CSV_PATH) || error("Missing $CSV_PATH — run paper_drf_validation.jl first.")

df  = CSV.read(CSV_PATH, DataFrame)
df.mass_kg = df.A_in2 .* df.L_in .* convert_to_m[:in]^3 .* ρ_STEEL

df_nlp = filter(r -> startswith(r.regime, "NLP"),       df)
df_cmp = filter(r -> startswith(r.regime, "Composite"), df)

@assert nrow(df_nlp) == nrow(df_cmp) "per-beam dataframe mismatch"
sort!(df_nlp, :beam_idx); sort!(df_cmp, :beam_idx)

m_nlp = df_nlp.mass_kg
m_cmp = df_cmp.mass_kg
n     = length(m_nlp)

println("Per-beam mass (n = $n beams, BaU r6c4):")
println("  total  — NLP: $(round(sum(m_nlp), digits=1)) kg, " *
        "Composite: $(round(sum(m_cmp), digits=1)) kg, " *
        "ratio: $(round(sum(m_nlp) / sum(m_cmp), digits=3))")
println("  per-beam m_nlp / m_cmp  summary:")
per_beam_ratio = m_nlp ./ m_cmp
for (k, v) in (("min", minimum(per_beam_ratio)),
               ("p25", quantile(per_beam_ratio, 0.25)),
               ("median", median(per_beam_ratio)),
               ("mean", mean(per_beam_ratio)),
               ("p75", quantile(per_beam_ratio, 0.75)),
               ("max", maximum(per_beam_ratio)))
    println("    $k  = $(round(v, digits=3))")
end

# ── figure ───────────────────────────────────────────────────────────────────
fontsize       = 14
fontsize_small = 11

fig = Figure(size = (460, 430), fontsize = fontsize)

axA = Axis(fig[1, 1],
    xlabel = "Nov 2024 NLP beam mass [kg]",
    ylabel = "Composite beam mass [kg]",
    title  = "Per-beam mass, BaU (r6c4, n = $n)",
    xgridvisible = false, ygridvisible = false,
    titlealign = :left,
    aspect = AxisAspect(1),
)

m_max = max(maximum(m_nlp), maximum(m_cmp)) * 1.05
lines!(axA, [0, m_max], [0, m_max];
       color = :black, linestyle = :dash, linewidth = 1, label = "y = x (parity)")

# colour perimeter beams differently — these carry half-tributary so usually
# they cluster near the lower-left of the scatter.
colour = [p ? 色[:lilac] : 色[:ceruleanblue] for p in df_nlp.is_perimeter]

scatter!(axA, m_nlp, m_cmp;
         color = colour, markersize = 7, alpha = 0.75)

xlims!(axA, 0, m_max); ylims!(axA, 0, m_max)

text!(axA, 0.04 * m_max, 0.94 * m_max;
      text  = "∑mₙₗₚ = $(round(sum(m_nlp), digits=0)) kg\n" *
              "∑m_cmp = $(round(sum(m_cmp), digits=0)) kg\n" *
              "ratio  = $(round(sum(m_nlp) / sum(m_cmp), digits=2))",
      align = (:left, :top), fontsize = fontsize_small)

# simple colour legend for interior / perimeter beams
elem_int   = MarkerElement(color = 色[:ceruleanblue], marker = :circle, markersize = 8)
elem_per   = MarkerElement(color = 色[:lilac],       marker = :circle, markersize = 8)
elem_par   = LineElement(color = :black, linestyle = :dash, linewidth = 1)
Legend(fig[2, 1], [elem_int, elem_per, elem_par],
       ["interior beam", "perimeter beam", "y = x (parity)"];
       orientation = :horizontal, tellheight = true, tellwidth = false,
       framevisible = false, labelsize = fontsize_small, patchsize = (14, 10))

save(FIG_PATH, fig)
println("\nWrote figure → $FIG_PATH")
