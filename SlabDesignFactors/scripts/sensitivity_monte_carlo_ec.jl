"""
Monte Carlo sensitivity study on embodied carbon coefficients.

Reads existing results CSVs (topology, grid, nova) and re-computes total EC
under 10,000 random draws of ECC_STEEL, ECC_CONCRETE, and ECC_REBAR from
log-normal distributions (CV ≈ 10%). No structural re-analysis is needed
because EC is a pure post-hoc calculation: `norm_mass * coefficient`.

Outputs:
  - Per-layout 90% confidence intervals on total EC
  - Ranking stability analysis (how often each layout is best/worst)
  - Tornado chart showing which coefficient dominates variance
  - CSV of all results

Addresses R2 #13, R4 #8 (uncertainty in EC coefficients).
"""

include("_scripts.jl")
using Random

CairoMakie.activate!()

# ── configuration ─────────────────────────────────────────────────────────────
Random.seed!(42)

N_SAMPLES = 10_000

deflection = "yes"
slabmin    = "yes"
results_base = "SlabDesignFactors/results/processed_$(deflection)deflection_$(slabmin)slabmin/"

save_path = "SlabDesignFactors/results/sensitivity_monte_carlo_ec/"
mkpath(save_path)

# ── nominal coefficients (from constants.jl, loaded by _scripts.jl) ───────────
# ECC_STEEL, ECC_CONCRETE, ECC_REBAR are already defined as global consts.
ECC_STEEL_NOM    = ECC_STEEL     # 1.22  kgCO₂e/kg
ECC_CONCRETE_NOM = ECC_CONCRETE  # 0.152 kgCO₂e/kg
ECC_REBAR_NOM    = ECC_REBAR     # 0.854 kgCO₂e/kg

# ── log-normal sampling with CV ≈ 10% (no Distributions.jl dependency) ────────
# If X ~ LogNormal(μ_ln, σ_ln), then E[X] = exp(μ_ln + σ_ln²/2)
# and CV = sqrt(exp(σ_ln²) - 1). For CV = 0.10, σ_ln ≈ 0.0998.
cv = 0.10
σ_ln = sqrt(log(1 + cv^2))

"""Sample N values from LogNormal with the given nominal mean and log-scale σ."""
function rand_lognormal(nominal::Real, σ::Real, N::Int)
    μ_ln = log(nominal) - σ^2 / 2
    return exp.(μ_ln .+ σ .* randn(N))
end

# Verify the sampler
test_samples = rand_lognormal(1.0, σ_ln, 100_000)
println("LogNormal sampler check (nominal=1.0, CV=$cv):")
println("  sample mean = $(round(mean(test_samples), digits=4))  (expect ≈ 1.0)")
println("  sample CV   = $(round(std(test_samples)/mean(test_samples), digits=4))  (expect ≈ $cv)")

# ── load existing data ────────────────────────────────────────────────────────
df_all = assemble_data(results_base)

# Filter to valid rows (non-zero area)
df = filter(row -> row.area > 0 && !isnan(row.steel_norm), df_all)
println("\nLoaded $(nrow(df)) valid layout configurations from $(results_base)")

# ── sample EC coefficients ────────────────────────────────────────────────────
ecc_steel_samples    = rand_lognormal(ECC_STEEL_NOM,    σ_ln, N_SAMPLES)
ecc_concrete_samples = rand_lognormal(ECC_CONCRETE_NOM, σ_ln, N_SAMPLES)
ecc_rebar_samples    = rand_lognormal(ECC_REBAR_NOM,    σ_ln, N_SAMPLES)

# ── Phase 2: recompute EC for every layout × every sample ─────────────────────
n_layouts = nrow(df)

# Matrix: rows = layouts, cols = MC samples
ec_steel_mat    = df.steel_norm    * ecc_steel_samples'
ec_concrete_mat = df.concrete_norm * ecc_concrete_samples'
ec_rebar_mat    = df.rebar_norm    * ecc_rebar_samples'
ec_total_mat    = ec_steel_mat .+ ec_concrete_mat .+ ec_rebar_mat

# ── per-layout statistics ─────────────────────────────────────────────────────
mc_results = DataFrame(
    name        = df.name,
    category    = df.category,
    rowcol      = df.rowcol,
    slab_type   = df.slab_type,
    beam_sizer  = df.beam_sizer,
    collinear   = df.collinear,
    steel_norm  = df.steel_norm,
    concrete_norm = df.concrete_norm,
    rebar_norm  = df.rebar_norm,
    ec_nominal  = df.total_ec,
    ec_mean     = vec(mean(ec_total_mat, dims=2)),
    ec_std      = vec(std(ec_total_mat, dims=2)),
    ec_p5       = [quantile(ec_total_mat[i, :], 0.05) for i in 1:n_layouts],
    ec_p25      = [quantile(ec_total_mat[i, :], 0.25) for i in 1:n_layouts],
    ec_p50      = [quantile(ec_total_mat[i, :], 0.50) for i in 1:n_layouts],
    ec_p75      = [quantile(ec_total_mat[i, :], 0.75) for i in 1:n_layouts],
    ec_p95      = [quantile(ec_total_mat[i, :], 0.95) for i in 1:n_layouts],
    ec_ci_width = [quantile(ec_total_mat[i, :], 0.95) - quantile(ec_total_mat[i, :], 0.05) for i in 1:n_layouts],
)

csv_path = joinpath(save_path, "monte_carlo_ec_results.csv")
CSV.write(csv_path, mc_results)
println("✓ Per-layout MC results saved to $csv_path ($(nrow(mc_results)) rows)")

# ── ranking stability ─────────────────────────────────────────────────────────
# For each MC sample, rank layouts by total EC (1 = lowest)
ranking_mat = zeros(Int, n_layouts, N_SAMPLES)
for j in 1:N_SAMPLES
    ranking_mat[:, j] = sortperm(sortperm(ec_total_mat[:, j]))
end

mc_results.rank_mean   = vec(mean(ranking_mat, dims=2))
mc_results.rank_median = [quantile(Float64.(ranking_mat[i, :]), 0.50) for i in 1:n_layouts]
mc_results.rank_std    = vec(std(ranking_mat, dims=2))
mc_results.pct_top5    = [sum(ranking_mat[i, :] .<= 5) / N_SAMPLES * 100 for i in 1:n_layouts]
mc_results.pct_bottom5 = [sum(ranking_mat[i, :] .>= n_layouts - 4) / N_SAMPLES * 100 for i in 1:n_layouts]

# Overwrite with ranking columns
CSV.write(csv_path, mc_results)
println("✓ Updated with ranking stability columns")

# ── relative reduction CI (vs BaU layout) ─────────────────────────────────────
# BaU = r1c2, uniaxial, discrete, collinear, uniform
bau_idx = findfirst(
    row -> row.name == "r1c2" &&
           row.slab_type == "uniaxial" &&
           row.beam_sizer == "discrete" &&
           row.collinear == true,
    eachrow(df)
)

if !isnothing(bau_idx)
    bau_ec_samples = ec_total_mat[bau_idx, :]
    println("\n── Reduction vs BaU ($(df.rowcol[bau_idx])) ──")
    println("  BaU nominal EC: $(round(df.total_ec[bau_idx], digits=2)) kgCO₂e/m²")
    println("  BaU MC:  $(round(mc_results.ec_p5[bau_idx], digits=2)) – $(round(mc_results.ec_p95[bau_idx], digits=2)) (90% CI)")

    # Best layout by nominal
    best_idx = argmin(df.total_ec)
    best_ec_samples = ec_total_mat[best_idx, :]
    reduction_samples = (bau_ec_samples .- best_ec_samples) ./ bau_ec_samples .* 100

    println("\n  Best layout: $(df.rowcol[best_idx])")
    println("  Nominal EC: $(round(df.total_ec[best_idx], digits=2)) kgCO₂e/m²")
    println("  Nominal reduction: $(round((df.total_ec[bau_idx] - df.total_ec[best_idx]) / df.total_ec[bau_idx] * 100, digits=1))%")
    println("  MC reduction: $(round(quantile(reduction_samples, 0.05), digits=1))% – $(round(quantile(reduction_samples, 0.95), digits=1))% (90% CI)")
    println("  MC mean reduction: $(round(mean(reduction_samples), digits=1))%")
end

# ── tornado / sensitivity decomposition ───────────────────────────────────────
# For each layout, compute variance contribution from each coefficient
# Using one-at-a-time approach: vary one, hold others at nominal
println("\n── Variance decomposition (% of total variance from each material) ──")

var_steel_vec    = Float64[]
var_concrete_vec = Float64[]
var_rebar_vec    = Float64[]

for i in 1:n_layouts
    s  = df.steel_norm[i]
    c  = df.concrete_norm[i]
    r  = df.rebar_norm[i]

    v_steel    = s^2 * var(ecc_steel_samples)
    v_concrete = c^2 * var(ecc_concrete_samples)
    v_rebar    = r^2 * var(ecc_rebar_samples)
    v_total    = v_steel + v_concrete + v_rebar

    if v_total > 0
        push!(var_steel_vec,    v_steel / v_total * 100)
        push!(var_concrete_vec, v_concrete / v_total * 100)
        push!(var_rebar_vec,    v_rebar / v_total * 100)
    else
        push!(var_steel_vec,    0.0)
        push!(var_concrete_vec, 0.0)
        push!(var_rebar_vec,    0.0)
    end
end

println("  Avg variance from steel:    $(round(mean(var_steel_vec), digits=1))%")
println("  Avg variance from concrete: $(round(mean(var_concrete_vec), digits=1))%")
println("  Avg variance from rebar:    $(round(mean(var_rebar_vec), digits=1))%")

# ── plots ─────────────────────────────────────────────────────────────────────
色 = Dict(
    :ceruleanblue => colorant"#00AEEF",
    :charcoalgrey => colorant"#3e3e3e",
    :magenta      => colorant"#dc267f",
    :gold         => colorant"#df7e00",
    :lilac        => colorant"#A678B5",
)

fontsize      = 11
smallfontsize = 8

# ── Plot 1: Error-bar chart (sorted by nominal EC) ───────────────────────────
sorted_idx = sortperm(mc_results.ec_nominal)
n_show = min(n_layouts, 60)
show_idx = sorted_idx[1:n_show]

fig1 = Figure(size=(1000, 500), fontsize=fontsize)
ax1 = Axis(fig1[1, 1],
    xlabel = "Layout (sorted by nominal EC)",
    ylabel = "Total EC (kgCO₂e/m²)",
    title  = "Monte Carlo 90% CI on embodied carbon (N=$N_SAMPLES, CV=$cv)",
    xticklabelsize = 6,
    xticklabelrotation = π/2,
    xticks = (1:n_show, mc_results.rowcol[show_idx]),
)

rangebars!(ax1, 1:n_show,
    mc_results.ec_p5[show_idx],
    mc_results.ec_p95[show_idx],
    color = (:grey, 0.4),
    linewidth = 1.5,
)
scatter!(ax1, 1:n_show, mc_results.ec_nominal[show_idx],
    color = 色[:ceruleanblue], markersize = 5, label = "Nominal",
)
scatter!(ax1, 1:n_show, mc_results.ec_mean[show_idx],
    color = 色[:magenta], markersize = 4, marker = :diamond, label = "MC mean",
)

axislegend(ax1, position = :lt, labelsize = smallfontsize, framevisible = false)
save(joinpath(save_path, "monte_carlo_ci.pdf"), fig1)
println("\n✓ CI plot saved")
display(fig1)

# ── Plot 2: Tornado chart (variance decomposition) ───────────────────────────
fig2 = Figure(size=(500, 400), fontsize=fontsize)
ax2 = Axis(fig2[1, 1],
    xlabel = "% of total EC variance",
    title  = "Variance decomposition by material",
    yticks = (1:3, ["Steel (ECC=$ECC_STEEL_NOM)", "Concrete (ECC=$ECC_CONCRETE_NOM)", "Rebar (ECC=$ECC_REBAR_NOM)"]),
)

means = [mean(var_steel_vec), mean(var_concrete_vec), mean(var_rebar_vec)]
colors = [色[:ceruleanblue], 色[:charcoalgrey], 色[:gold]]
barplot!(ax2, 1:3, means, direction = :x, color = colors, bar_labels = :values,
    label_formatter = x -> "$(round(x, digits=1))%")
xlims!(ax2, (0, 100))

save(joinpath(save_path, "monte_carlo_tornado.pdf"), fig2)
println("✓ Tornado plot saved")
display(fig2)

# ── Plot 3: Ranking stability histogram for top layouts ──────────────────────
fig3 = Figure(size=(700, 400), fontsize=fontsize)
ax3 = Axis(fig3[1, 1],
    xlabel = "Rank (1 = lowest EC)",
    ylabel = "Frequency",
    title  = "Rank distribution for top-5 nominal layouts (N=$N_SAMPLES)",
)

nominal_top5 = sorted_idx[1:min(5, n_layouts)]
for (k, i) in enumerate(nominal_top5)
    hist!(ax3, Float64.(ranking_mat[i, :]),
        bins = 30,
        label = mc_results.rowcol[i],
        color = (colors[mod1(k, length(colors))], 0.5),
        strokewidth = 0.5,
    )
end

axislegend(ax3, position = :rt, labelsize = smallfontsize, framevisible = false)
save(joinpath(save_path, "monte_carlo_rank_stability.pdf"), fig3)
println("✓ Ranking stability plot saved")
display(fig3)

# ── Plot 4: Reduction CI vs BaU ──────────────────────────────────────────────
if !isnothing(bau_idx)
    fig4 = Figure(size=(700, 400), fontsize=fontsize)
    ax4 = Axis(fig4[1, 1],
        xlabel = "EC reduction vs BaU (%)",
        ylabel = "Frequency",
        title  = "Distribution of EC reduction for best layout vs BaU ($(df.rowcol[bau_idx]))",
    )

    best_idx_local = argmin(df.total_ec)
    reduction_samples_plot = (bau_ec_samples .- ec_total_mat[best_idx_local, :]) ./ bau_ec_samples .* 100

    hist!(ax4, reduction_samples_plot,
        bins = 50,
        color = (色[:ceruleanblue], 0.7),
        strokewidth = 0.5,
    )
    vlines!(ax4, [mean(reduction_samples_plot)], color = 色[:magenta], linewidth = 2, label = "Mean")
    vlines!(ax4, [quantile(reduction_samples_plot, 0.05), quantile(reduction_samples_plot, 0.95)],
        color = 色[:gold], linewidth = 1.5, linestyle = :dash, label = "90% CI",
    )

    axislegend(ax4, position = :lt, labelsize = smallfontsize, framevisible = false)
    save(joinpath(save_path, "monte_carlo_reduction_ci.pdf"), fig4)
    println("✓ Reduction CI plot saved")
    display(fig4)
end

println("\n════════════════════════════════════════════")
println("  Monte Carlo EC sensitivity complete")
println("  $N_SAMPLES samples, CV = $(cv * 100)%")
println("  Results: $save_path")
println("════════════════════════════════════════════")
