"""
Monte Carlo sensitivity study on embodied carbon coefficients, over three
deterministic steel scenarios (low / typical / high steel EC intensity).

Motivation
----------
To address the "uniform material properties" reviewer comment, we bracket
steel-type uncertainty with three center points (kgCO₂e/kg, cradle-to-gate):

    Low:     CLF 2021 P20 for hot-rolled structural sections (Carlisle et al. 2021).
    Typical: CLF 2023 fabricated hot-rolled sections baseline (Waldman et al. 2023,
             AISC 2021 IW-EPD) — aligned with `ECC_STEEL = 1.22` in the main pipeline.
    High:    CLF 2021 P80 for hot-rolled structural sections (Carlisle et al. 2021).

Low and high use CLF 2021 percentiles as plausible market bounds; CLF 2023
Appendix D4 does not publish percentiles (insufficient product EPDs). The MC
study remains self-consistent: at each scenario center, every sample uses that
center for the steel coefficient mean before log-normal dispersion.

At each steel center, we run 10,000 log-normal draws (CV = 10%) on all three
coefficients (steel, concrete, rebar). Total EC accounting is steel +
concrete + rebar only (no fireproofing or columns), matching the underlying
processed CSVs. No structural re-analysis is required: EC is a pure post-hoc
linear calculation `norm_mass × coefficient`.

The sweep is additionally run over two *design cases* — `nodeflection_yesslabmin`
(no deflection check, minimum slab thickness enforced, BaU-aligned) and
`nodeflection_noslabmin` (no deflection check, slab-min relaxed) — since
relaxing the minimum-slab constraint changes the best-layout identity and
therefore the headline reduction number.

Outputs per design case (subfolder `<design_case>/`), per steel scenario
(sub-subfolder `low/`, `typical/`, `high/`):
  - `monte_carlo_ec_results.csv`     per-layout stats (nominal, mean, std, P5/P25/P50/P75/P95, ranking)
  - `monte_carlo_caterpillar.pdf`    EC distribution across all layouts (shaded 90% band, best + BaU annotated)
  - `monte_carlo_tornado.pdf`        variance decomposition by material
  - `monte_carlo_reduction_ci.pdf`   reduction-vs-BaU distribution with mean + 90% CI

Per design case (at `save_path_root/<design_case>/`):
  - `scenario_summary.csv`                          per-scenario BaU/best/reduction stats (paper-ready table)
  - `forest_reduction_by_steel_scenario.pdf`        HEADLINE: reduction ± 90% CI by scenario with best-layout identity
  - `violin_reduction_by_steel_scenario.pdf`        reduction distributions side-by-side
  - `rank_stability_summary.pdf`                    % of MC samples in which nominal-best stays rank 1

Cross-design outputs (at `save_path_root/`):
  - `cross_design_summary.csv`                      all (design_case × steel_scenario) rows in one table

References:
- Waldman, B., Hyatt, A., Carlisle, S., Palmeri, J., & Simonen, K. (2023).
  *2023 CLF North American Material Baselines* (Structural Steel, Appendix D4) — Typical center.
- Carlisle, S. et al. (2021). *2021 CLF Material Baseline Report*
  (Hot-Rolled Structural Steel Sections) — P20/P80 used for Low/High bounds.
- Hart, J., D'Amico, B., & Pomponi, F. (2021). *Whole-life embodied carbon in multistory
  buildings*, J. Ind. Ecol. — precedent for log-normal ECC distributions (CV ≈ 10%).

Addresses R2 #13, R4 #8 (uncertainty in EC coefficients) and the
"uniform material properties" reviewer comment (steel-type variability).
"""

include("_scripts.jl")
using Random

CairoMakie.activate!()

# ── configuration ─────────────────────────────────────────────────────────────
N_SAMPLES = 10_000
SEED      = 42

# Design cases — each pairs a deflection flag with a slab-min flag and points
# at its own processed-results folder. Running both surfaces how relaxing the
# minimum-slab constraint shifts the best-layout identity and the headline
# reduction distribution.
design_cases = [
    (name         = "nodeflection_yesslabmin",
     label        = "Slab minimum of 0.125m",
     results_base = "SlabDesignFactors/results/processed_nodeflection_yesslabmin/"),
    (name         = "nodeflection_noslabmin",
     label        = "No slab minimum",
     results_base = "SlabDesignFactors/results/processed_nodeflection_noslabmin/"),
]

save_path_root = "SlabDesignFactors/results/sensitivity_monte_carlo_ec/"
mkpath(save_path_root)

# Log-normal sampling with CV ≈ 10%: σ_ln = sqrt(log(1 + CV²)) ≈ 0.0998
const CV   = 0.10
const σ_ln = sqrt(log(1 + CV^2))

# Steel scenario center points (kgCO₂e/kg).
# Typical: CLF 2023 fabricated hot-rolled sections (Waldman et al. 2023, AISC 2021 IW-EPD).
# Low/High: CLF 2021 P20/P80 for hot-rolled structural sections (Carlisle et al. 2021),
#           used as plausible market bounds since CLF 2023 Appendix D4 does not publish
#           percentiles (insufficient product EPDs).
steel_scenarios = [
    (name = "low",     ecc_steel = 0.80, label = "Low (CLF 2021 P20)"),
    (name = "typical", ecc_steel = 1.22, label = "Typical (CLF 2023 baseline)"),
    (name = "high",    ecc_steel = 1.70, label = "High (CLF 2021 P80)"),
]

# ── helpers ───────────────────────────────────────────────────────────────────
"""Sample N values from LogNormal whose mean equals `nominal` (CV controlled by σ)."""
function rand_lognormal(nominal::Real, σ::Real, N::Int)
    μ_ln = log(nominal) - σ^2 / 2
    return exp.(μ_ln .+ σ .* randn(N))
end

"""
    find_bau_idx(df) -> Union{Int, Nothing}

Locate the canonical business-as-usual (BaU) slab layout using the *full*
paper filter (matching `plot/plotting/8_stats_summary.jl`,
`plot/plotting/9_stats_topology.jl`, and `plot/plotting/2_megaplot.jl`):

    name == "r1c2" ∧ slab_type == "uniaxial" ∧ slab_sizer == "uniform"
    ∧ beam_sizer == "discrete" ∧ collinear == true
    ∧ vector_1d_x == 1 ∧ vector_1d_y == 0 ∧ max_depth == 40

Earlier versions used only the first four criteria, which matched multiple
rows; `findfirst` then picked whichever `r1c2` variant appeared first after
`sort!(df, [:row, :col])`, giving a BaU EC ~1.2% lower than the canonical
value and a correspondingly inflated reduction.
"""
function find_bau_idx(df)
    return findfirst(
        row -> row.name         == "r1c2"     &&
               row.slab_type    == "uniaxial" &&
               row.slab_sizer   == "uniform"  &&
               row.beam_sizer   == "discrete" &&
               row.collinear    == true       &&
               row.vector_1d_x  == 1          &&
               row.vector_1d_y  == 0          &&
               row.max_depth    == 40,
        eachrow(df),
    )
end

"""
    run_mc_scenario(df, ecc_steel_center; N, σ, seed) -> (results, ec_total_mat, ranking_mat, var_contrib)

Runs one scenario (fixed steel center, log-normal MC on all three coefficients)
and returns a per-layout statistics DataFrame, the full EC sample matrix
(layouts × samples), the ranking matrix, and the average per-material variance
contribution (%) across layouts.
"""
function run_mc_scenario(df, ecc_steel_center::Real;
                         N::Int = N_SAMPLES, σ::Real = σ_ln, seed::Int = SEED)
    Random.seed!(seed)

    ecc_steel    = rand_lognormal(ecc_steel_center, σ, N)
    ecc_concrete = rand_lognormal(ECC_CONCRETE,     σ, N)
    ecc_rebar    = rand_lognormal(ECC_REBAR,        σ, N)

    n = nrow(df)

    ec_total_mat = df.steel_norm    * ecc_steel'    .+
                   df.concrete_norm * ecc_concrete' .+
                   df.rebar_norm    * ecc_rebar'

    ec_nominal = df.steel_norm    .* ecc_steel_center .+
                 df.concrete_norm .* ECC_CONCRETE     .+
                 df.rebar_norm    .* ECC_REBAR

    ranking_mat = zeros(Int, n, N)
    for j in 1:N
        ranking_mat[:, j] = sortperm(sortperm(ec_total_mat[:, j]))
    end

    results = DataFrame(
        name          = df.name,
        category      = df.category,
        rowcol        = df.rowcol,
        slab_type     = df.slab_type,
        beam_sizer    = df.beam_sizer,
        collinear     = df.collinear,
        steel_norm    = df.steel_norm,
        concrete_norm = df.concrete_norm,
        rebar_norm    = df.rebar_norm,
        ec_nominal    = ec_nominal,
        ec_mean       = vec(mean(ec_total_mat, dims = 2)),
        ec_std        = vec(std(ec_total_mat,  dims = 2)),
        ec_p5         = [quantile(ec_total_mat[i, :], 0.05) for i in 1:n],
        ec_p25        = [quantile(ec_total_mat[i, :], 0.25) for i in 1:n],
        ec_p50        = [quantile(ec_total_mat[i, :], 0.50) for i in 1:n],
        ec_p75        = [quantile(ec_total_mat[i, :], 0.75) for i in 1:n],
        ec_p95        = [quantile(ec_total_mat[i, :], 0.95) for i in 1:n],
        ec_ci_width   = [quantile(ec_total_mat[i, :], 0.95) -
                         quantile(ec_total_mat[i, :], 0.05) for i in 1:n],
        rank_mean     = vec(mean(ranking_mat, dims = 2)),
        rank_median   = [quantile(Float64.(ranking_mat[i, :]), 0.50) for i in 1:n],
        rank_std      = vec(std(ranking_mat, dims = 2)),
        pct_top5      = [sum(ranking_mat[i, :] .<= 5) / N * 100 for i in 1:n],
        pct_bottom5   = [sum(ranking_mat[i, :] .>= n - 4) / N * 100 for i in 1:n],
    )

    # Variance decomposition (first-order, one-at-a-time), averaged across layouts
    v_steel    = df.steel_norm    .^ 2 .* var(ecc_steel)
    v_concrete = df.concrete_norm .^ 2 .* var(ecc_concrete)
    v_rebar    = df.rebar_norm    .^ 2 .* var(ecc_rebar)
    v_total    = v_steel .+ v_concrete .+ v_rebar
    safe_total = ifelse.(v_total .> 0, v_total, 1.0)
    var_contrib = (
        steel    = mean(v_steel    ./ safe_total) * 100,
        concrete = mean(v_concrete ./ safe_total) * 100,
        rebar    = mean(v_rebar    ./ safe_total) * 100,
    )

    return results, ec_total_mat, ranking_mat, var_contrib
end

# ── plotting helpers ──────────────────────────────────────────────────────────
# `色` palette is defined globally in `SlabDesignFactors/plot/plotting/utils.jl`
# and loaded via `_scripts.jl`; do not redeclare here (Julia 1.12 forbids
# promoting an existing global to const).

const FONTSIZE       = 14
const SMALLFONTSIZE  = 11
const SCENARIO_PALETTE = [色[:ceruleanblue], 色[:charcoalgrey], 色[:magenta]]

"""
    plot_scenario_caterpillar(df, results, label, outpath)

Caterpillar plot: every layout sorted by nominal EC, with the P5–P95 MC band
shaded around the nominal line. Best and BaU are annotated as points; no
per-layout x-labels (the distribution shape carries the message).
"""
function plot_scenario_caterpillar(df, results, label, outpath)
    n          = nrow(results)
    sorted_idx = sortperm(results.ec_nominal)
    xs         = 1:n
    p5         = results.ec_p5[sorted_idx]
    p95        = results.ec_p95[sorted_idx]
    nominal    = results.ec_nominal[sorted_idx]

    fig = Figure(size = (800, 450), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "Layouts sorted by nominal EC  (n = $n)",
        ylabel = "Total EC (kgCO₂e/m²)",
        title  = "$label  —  MC 90% band (N=$N_SAMPLES, CV=$CV)",
        xticksvisible      = false,
        xticklabelsvisible = false,
    )

    band!(ax, xs, p5, p95, color = (色[:ceruleanblue], 0.25),
        label = "90% CI")
    lines!(ax, xs, nominal, color = 色[:ceruleanblue], linewidth = 1.5,
        label = "Nominal")

    scatter!(ax, [1], [nominal[1]], color = 色[:gold], markersize = 16)
    text!(ax, "Best · $(results.rowcol[sorted_idx[1]])",
        position = (1, nominal[1]),
        offset = (12, -4), align = (:left, :center),
        fontsize = SMALLFONTSIZE)

    bau_idx = find_bau_idx(df)
    if !isnothing(bau_idx)
        bau_pos = findfirst(==(bau_idx), sorted_idx)
        if !isnothing(bau_pos)
            scatter!(ax, [bau_pos], [nominal[bau_pos]],
                color = 色[:magenta], markersize = 16)
            text!(ax, "BaU · $(df.rowcol[bau_idx])",
                position = (bau_pos, nominal[bau_pos]),
                offset = (12, 6), align = (:left, :center),
                fontsize = SMALLFONTSIZE)
        end
    end

    axislegend(ax, position = :lt, labelsize = SMALLFONTSIZE, framevisible = false)
    save(outpath, fig)
    return fig
end

function plot_scenario_tornado(var_contrib, label, ecc_steel_center, outpath)
    fig = Figure(size = (600, 400), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "% of total EC variance",
        title  = "Variance decomposition — $label",
        yticks = (1:3, [
            "Steel (ECC=$(round(ecc_steel_center, digits=2)))",
            "Concrete (ECC=$ECC_CONCRETE)",
            "Rebar (ECC=$ECC_REBAR)",
        ]),
    )
    values     = [var_contrib.steel, var_contrib.concrete, var_contrib.rebar]
    bar_colors = [色[:ceruleanblue], 色[:charcoalgrey], 色[:gold]]
    bar_labels = ["$(round(v, digits=1))%" for v in values]
    barplot!(ax, 1:3, values, direction = :x, color = bar_colors,
        bar_labels = bar_labels)
    xlims!(ax, (0, 110))
    save(outpath, fig)
    return fig
end

function plot_scenario_reduction_ci(df, ec_total_mat, results, label, outpath)
    bau_idx = find_bau_idx(df)
    isnothing(bau_idx) && return nothing

    best_idx  = argmin(results.ec_nominal)
    bau_samp  = ec_total_mat[bau_idx,  :]
    best_samp = ec_total_mat[best_idx, :]
    reduction = (bau_samp .- best_samp) ./ bau_samp .* 100

    fig = Figure(size = (700, 400), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "EC reduction vs BaU (%)",
        ylabel = "Frequency",
        title  = "Reduction best-vs-BaU ($(df.rowcol[bau_idx])) — $label",
    )
    hist!(ax, reduction, bins = 50, color = (色[:ceruleanblue], 0.7), strokewidth = 0.5)
    vlines!(ax, [mean(reduction)],
        color = 色[:magenta], linewidth = 2, label = "Mean")
    vlines!(ax, [quantile(reduction, 0.05), quantile(reduction, 0.95)],
        color = 色[:gold], linewidth = 1.5, linestyle = :dash, label = "90% CI")
    axislegend(ax, position = :lt, labelsize = SMALLFONTSIZE, framevisible = false)
    save(outpath, fig)
    return reduction
end

"""
    plot_forest(summary_records, steel_scenarios, outpath)

Headline MC figure: horizontal forest plot with nominal reduction ± 90% CI per
steel scenario. Labels are sourced from `steel_scenarios` (keyed by scenario
name) and embed the best-layout identity so ranking stability is visible at a
glance.
"""
function plot_forest(summary_records, steel_scenarios, outpath)
    n = nrow(summary_records)
    y = collect(n:-1:1)  # top-to-bottom reading order follows scenario list

    label_map = Dict(s.name => s.label for s in steel_scenarios)
    y_labels = ["$(label_map[row.scenario])\nECC_STEEL=$(row.ecc_steel_center)\nbest=$(row.best_rowcol)"
                for row in eachrow(summary_records)]

    fig = Figure(size = (900, 450), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "EC reduction best-vs-BaU (%)",
        yticks = (y, y_labels),
        title  = "MC embodied-carbon savings across steel scenarios  (N=$N_SAMPLES, CV=$CV)",
        yticklabelsize = SMALLFONTSIZE,
    )

    for (k, row) in enumerate(eachrow(summary_records))
        rangebars!(ax, [y[k]], [row.reduction_p5], [row.reduction_p95],
            direction = :x, color = SCENARIO_PALETTE[k],
            linewidth = 3, whiskerwidth = 14)
        scatter!(ax, [row.reduction_nominal], [y[k]],
            color = SCENARIO_PALETTE[k], markersize = 16, strokewidth = 1, strokecolor = :black)
        text!(ax,
            "$(round(row.reduction_nominal, digits=1))%  [$(round(row.reduction_p5, digits=1))–$(round(row.reduction_p95, digits=1))]",
            position = (row.reduction_p95, y[k]),
            offset = (10, 0), align = (:left, :center),
            fontsize = SMALLFONTSIZE)
    end

    # Pad x so the annotation text doesn't clip
    xmin = minimum(summary_records.reduction_p5)  - 2
    xmax = maximum(summary_records.reduction_p95) + 14
    xlims!(ax, (xmin, xmax))
    ylims!(ax, (0.4, n + 0.6))
    save(outpath, fig)
    return fig
end

"""
    dump_reduction_csvs(df, scenario_store, steel_scenarios, save_path)

Persist the full best-vs-BaU reduction samples and per-scenario summary stats
to CSV so the violin plot can be iterated on via `replot_monte_carlo_violin.jl`
without rerunning the Monte Carlo loop.

Writes two files into `save_path`:
  • `reduction_samples.csv`   wide format, one column per scenario name, N rows
  • `reduction_summary.csv`   per-scenario nominal / mean / P5-P25-P50-P75-P95
                              plus labels and the BaU / best layout identities
"""
function dump_reduction_csvs(df, scenario_store, steel_scenarios, save_path)
    bau_idx = find_bau_idx(df)
    isnothing(bau_idx) && return nothing

    samples_df   = DataFrame()
    summary_rows = NamedTuple[]

    for scenario in steel_scenarios
        r         = scenario_store[scenario.name]
        best_idx  = argmin(r.results.ec_nominal)
        bau_samp  = r.ec_total_mat[bau_idx,  :]
        best_samp = r.ec_total_mat[best_idx, :]
        reduction = (bau_samp .- best_samp) ./ bau_samp .* 100

        samples_df[!, scenario.name] = reduction

        p5, p25, p50, p75, p95 = quantile(reduction, (0.05, 0.25, 0.50, 0.75, 0.95))
        nominal = (r.results.ec_nominal[bau_idx] - r.results.ec_nominal[best_idx]) /
                   r.results.ec_nominal[bau_idx] * 100

        push!(summary_rows, (
            scenario         = scenario.name,
            label            = scenario.label,
            ecc_steel_center = scenario.ecc_steel,
            bau_rowcol       = df.rowcol[bau_idx],
            best_rowcol      = r.results.rowcol[best_idx],
            nominal          = nominal,
            mean             = mean(reduction),
            std              = std(reduction),
            p5               = p5,
            p25              = p25,
            p50              = p50,
            p75              = p75,
            p95              = p95,
        ))
    end

    summary_df = DataFrame(summary_rows)

    samples_path = joinpath(save_path, "reduction_samples.csv")
    summary_path = joinpath(save_path, "reduction_summary.csv")
    CSV.write(samples_path, samples_df)
    CSV.write(summary_path, summary_df)

    return samples_path, summary_path
end

"""
    plot_rank_stability_summary(scenario_store, steel_scenarios, outpath)

Single bar chart answering: "in what % of MC samples does the nominal-best
layout remain rank 1?" — one bar per scenario.
"""
function plot_rank_stability_summary(scenario_store, steel_scenarios, outpath)
    bars   = Float64[]
    labels = String[]
    for scenario in steel_scenarios
        r          = scenario_store[scenario.name]
        best_idx   = argmin(r.results.ec_nominal)
        pct_rank_1 = sum(r.ranking_mat[best_idx, :] .== 1) / size(r.ranking_mat, 2) * 100
        push!(bars, pct_rank_1)
        push!(labels, "$(scenario.label)\nbest=$(r.results.rowcol[best_idx])")
    end

    fig = Figure(size = (700, 420), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        ylabel = "% of MC samples with nominal-best ranked #1",
        xticks = (1:length(bars), labels),
        xticklabelsize = SMALLFONTSIZE,
        title = "Ranking stability of the nominal-best layout",
    )
    barplot!(ax, 1:length(bars), bars,
        color = SCENARIO_PALETTE[1:length(bars)],
        bar_labels = ["$(round(b, digits=1))%" for b in bars])
    ylims!(ax, (0, 110))
    save(outpath, fig)
    return fig
end

# Include the violin replot helper once up front; re-including inside the loop
# would just be wasteful I/O.
include("replot_monte_carlo_violin.jl")

# Sampler sanity check (once — independent of design case)
test_samples = rand_lognormal(1.0, σ_ln, 100_000)
println("LogNormal sampler check (nominal=1.0, CV=$CV):")
println("  sample mean = $(round(mean(test_samples), digits=4))  (expect ≈ 1.0)")
println("  sample CV   = $(round(std(test_samples)/mean(test_samples), digits=4))  (expect ≈ $CV)")

"""
    run_design_case(design_case, save_path_root) -> DataFrame

Run the full Monte Carlo × steel-scenario sweep for a single design case
(deflection × slab-min combination). Writes per-scenario and cross-scenario
outputs to `save_path_root/<design_case.name>/`. Returns the
`scenario_summary` DataFrame tagged with a `design_case` column so the caller
can concatenate across design cases.
"""
function run_design_case(design_case, save_path_root)
    println("\n████████████████████████████████████████████")
    println("  Design case: $(design_case.label)")
    println("  Data:        $(design_case.results_base)")
    println("████████████████████████████████████████████")

    save_path = joinpath(save_path_root, design_case.name)
    mkpath(save_path)

    # ── load data ────────────────────────────────────────────────────────────
    df_all = assemble_data(design_case.results_base)
    df     = filter(row -> row.area > 0 && !isnan(row.steel_norm), df_all)
    println("\nLoaded $(nrow(df)) valid layout configurations from $(design_case.results_base)")

    # ── per-scenario loop ────────────────────────────────────────────────────
    scenario_store  = Dict{String, Any}()
    summary_records = DataFrame(
        scenario           = String[],
        ecc_steel_center   = Float64[],
        bau_rowcol         = String[],
        bau_ec_nominal     = Float64[],
        bau_ec_p5          = Float64[],
        bau_ec_p95         = Float64[],
        best_rowcol        = String[],
        best_ec_nominal    = Float64[],
        best_ec_p5         = Float64[],
        best_ec_p95        = Float64[],
        reduction_nominal  = Float64[],
        reduction_mean     = Float64[],
        reduction_p5       = Float64[],
        reduction_p95      = Float64[],
        var_pct_steel      = Float64[],
        var_pct_concrete   = Float64[],
        var_pct_rebar      = Float64[],
    )

    for scenario in steel_scenarios
        println("\n════════════════════════════════════════════")
        println("  Steel scenario: $(scenario.label)   ECC_STEEL = $(scenario.ecc_steel)")
        println("════════════════════════════════════════════")

        scenario_dir = joinpath(save_path, scenario.name)
        mkpath(scenario_dir)

        results, ec_total_mat, ranking_mat, var_contrib =
            run_mc_scenario(df, scenario.ecc_steel)

        csv_path = joinpath(scenario_dir, "monte_carlo_ec_results.csv")
        CSV.write(csv_path, results)
        println("  ✓ Results saved → $csv_path  ($(nrow(results)) rows)")

        plot_scenario_caterpillar(df, results, scenario.label,
            joinpath(scenario_dir, "monte_carlo_caterpillar.pdf"))
        plot_scenario_tornado(var_contrib, scenario.label, scenario.ecc_steel,
            joinpath(scenario_dir, "monte_carlo_tornado.pdf"))
        plot_scenario_reduction_ci(df, ec_total_mat, results, scenario.label,
            joinpath(scenario_dir, "monte_carlo_reduction_ci.pdf"))

        bau_idx = find_bau_idx(df)
        if !isnothing(bau_idx)
            best_idx  = argmin(results.ec_nominal)
            bau_samp  = ec_total_mat[bau_idx,  :]
            best_samp = ec_total_mat[best_idx, :]
            reduction = (bau_samp .- best_samp) ./ bau_samp .* 100

            push!(summary_records, (
                scenario          = scenario.name,
                ecc_steel_center  = scenario.ecc_steel,
                bau_rowcol        = df.rowcol[bau_idx],
                bau_ec_nominal    = results.ec_nominal[bau_idx],
                bau_ec_p5         = results.ec_p5[bau_idx],
                bau_ec_p95        = results.ec_p95[bau_idx],
                best_rowcol       = df.rowcol[best_idx],
                best_ec_nominal   = results.ec_nominal[best_idx],
                best_ec_p5        = results.ec_p5[best_idx],
                best_ec_p95       = results.ec_p95[best_idx],
                reduction_nominal = (results.ec_nominal[bau_idx] - results.ec_nominal[best_idx]) /
                                     results.ec_nominal[bau_idx] * 100,
                reduction_mean    = mean(reduction),
                reduction_p5      = quantile(reduction, 0.05),
                reduction_p95     = quantile(reduction, 0.95),
                var_pct_steel     = var_contrib.steel,
                var_pct_concrete  = var_contrib.concrete,
                var_pct_rebar     = var_contrib.rebar,
            ))

            println("  BaU $(df.rowcol[bau_idx])  nominal EC = $(round(results.ec_nominal[bau_idx], digits=2))")
            println("  Best $(df.rowcol[best_idx]) nominal EC = $(round(results.ec_nominal[best_idx], digits=2))")
            println("  Nominal reduction: $(round((results.ec_nominal[bau_idx] - results.ec_nominal[best_idx]) / results.ec_nominal[bau_idx] * 100, digits=1))%")
            println("  MC reduction:      $(round(quantile(reduction, 0.05), digits=1))% – $(round(quantile(reduction, 0.95), digits=1))% (90% CI)")
        end
        println("  Variance: steel $(round(var_contrib.steel, digits=1))%, concrete $(round(var_contrib.concrete, digits=1))%, rebar $(round(var_contrib.rebar, digits=1))%")

        scenario_store[scenario.name] = (; results, ec_total_mat, ranking_mat, var_contrib, scenario)
    end

    summary_csv = joinpath(save_path, "scenario_summary.csv")
    CSV.write(summary_csv, summary_records)
    println("\n✓ Cross-scenario summary → $summary_csv")

    # ── cross-scenario comparison plots ──────────────────────────────────────
    println("\n════════════════════════════════════════════")
    println("  Cross-scenario comparison  ($(design_case.label))")
    println("════════════════════════════════════════════")

    plot_forest(summary_records, steel_scenarios,
        joinpath(save_path, "forest_reduction_by_steel_scenario.pdf"))
    println("  ✓ Headline forest plot saved")

    samples_path, summary_path = dump_reduction_csvs(df, scenario_store, steel_scenarios, save_path)
    println("  ✓ Reduction samples → $samples_path")
    println("  ✓ Reduction summary → $summary_path")

    replot_monte_carlo_violin(save_path;
        outpath = joinpath(save_path, "violin_reduction_by_steel_scenario.pdf"))
    println("  ✓ Cross-scenario reduction violin saved")

    replot_monte_carlo_layoutspace(save_path;
        outpath = joinpath(save_path, "layoutspace_reduction_by_steel_scenario.pdf"))
    println("  ✓ Layout-space reduction plot saved (design variance + MC CI)")

    plot_rank_stability_summary(scenario_store, steel_scenarios,
        joinpath(save_path, "rank_stability_summary.pdf"))
    println("  ✓ Rank-stability summary saved")

    println("\nBest layout (lowest nominal EC) by scenario:")
    for scenario in steel_scenarios
        r        = scenario_store[scenario.name]
        best_idx = argmin(r.results.ec_nominal)
        println("  $(rpad(scenario.label, 20)) → $(r.results.rowcol[best_idx])   EC = $(round(r.results.ec_nominal[best_idx], digits=2))")
    end

    summary_records[!, :design_case] .= design_case.name
    return summary_records
end

# ── run each design case and aggregate ────────────────────────────────────────
all_summaries = DataFrame[]
for design_case in design_cases
    push!(all_summaries, run_design_case(design_case, save_path_root))
end

cross_design = reduce(vcat, all_summaries)
# Promote `design_case` to leftmost column for readability.
select!(cross_design, :design_case, Not(:design_case))

cross_csv = joinpath(save_path_root, "cross_design_summary.csv")
CSV.write(cross_csv, cross_design)
println("\n✓ Cross-design summary → $cross_csv")

# Combined violin: paired violins per steel scenario overlaying the design
# cases, so the reader sees "with vs without slab-min" at a glance.
combined_violin_path = joinpath(save_path_root, "violin_reduction_combined.pdf")
replot_monte_carlo_violin_combined(save_path_root, design_cases;
    outpath = combined_violin_path)
println("✓ Combined design-case violin → $combined_violin_path")

combined_layout_path = joinpath(save_path_root, "layoutspace_reduction_combined.pdf")
replot_monte_carlo_layoutspace_combined(save_path_root, design_cases;
    outpath = combined_layout_path)
println("✓ Combined design-case layout-space plot → $combined_layout_path")

# ── cross-design headline comparison ──────────────────────────────────────────
println("\n════════════════════════════════════════════")
println("  Cross-design comparison: best-layout identity")
println("════════════════════════════════════════════")
for design_case in design_cases
    rows = filter(row -> row.design_case == design_case.name, cross_design)
    println("\n  $(design_case.label)  ($(design_case.name))")
    for row in eachrow(rows)
        println("    $(rpad(row.scenario, 8)) (ECC_STEEL=$(row.ecc_steel_center))" *
                "  best=$(row.best_rowcol)  reduction=$(round(row.reduction_nominal, digits=1))%" *
                "  [90% CI $(round(row.reduction_p5, digits=1))–$(round(row.reduction_p95, digits=1))%]")
    end
end

println("\n════════════════════════════════════════════")
println("  Monte Carlo EC sensitivity complete")
println("  $(length(design_cases)) design cases × $(length(steel_scenarios)) steel scenarios × $N_SAMPLES samples (CV = $(CV*100)%)")
println("  Results: $save_path_root")
println("════════════════════════════════════════════")
