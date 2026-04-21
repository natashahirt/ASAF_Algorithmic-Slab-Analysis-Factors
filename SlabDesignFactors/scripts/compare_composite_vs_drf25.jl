"""
Compare composite-beam pipeline (composite_action=true) against two bare-steel
baselines on a representative subset of topology geometries:

  1. `composite` — composite_action=true, drf=1.0 (new staged pipeline).
  2. `naive25`   — bare steel, deflection constraints OFF during sizing
                   (strength-only optimization). Serviceability is then judged
                   only in post-processing, dividing computed δ by 2.5 as a
                   naive proxy for composite stiffening.

A legacy CSV baseline (remote_results_yesdeflection_yesslabmin/topology.csv),
produced by the dev-remote pipeline (bare steel with strict L/360 on the
bare section, i.e. effective drf=1.0), is also loaded for reference.

Usage (from project root):
    julia --project=. SlabDesignFactors/scripts/compare_composite_vs_drf25.jl
"""

include("_scripts.jl")
using PrettyTables

println("\n══════════════════════════════════════════════")
println("  Composite vs Naive-2.5× vs Legacy-L/360 Comparison Sweep")
println("══════════════════════════════════════════════\n")

# ── Configuration ─────────────────────────────────────────────────────────────

main_path = "Geometries/topology/"
old_csv_path = "SlabDesignFactors/results/remote_results_yesdeflection_yesslabmin/topology.csv"
results_path = "SlabDesignFactors/results/comparison_composite/"
results_name = "topology_composite"

geom_files = [
    "r1c1.json", "r1c2.json",
    "r2c3.json", "r3c2.json",
    "r5c2.json", "r6c4.json",
    "r7c4.json", "r9c4.json",
]

max_depths = [25, 40]
slab_sizer = :uniform
slab_type = :isotropic
vector_1d = [0.0, 0.0]

# Sizing variants evaluated per geometry/depth.
# - `composite` sizes with staged composite deflection (drf=1.0, strict L/360).
# - `naive25`   sizes bare steel with deflection constraints OFF (strength-only
#   optimization, drf=1.0 in the sizer is inert). Serviceability is then
#   judged only in post-processing via `apply_naive_drf`, which divides the
#   computed bare-beam δ by 2.5 as a naive composite-stiffening proxy.
variants = [
    (name = "composite", composite = true,  drf = 1.0, deflection_limit = true),
    (name = "naive25",   composite = false, drf = 1.0, deflection_limit = false),
]
const NAIVE_DRF = 2.5  # post-processing-only stiffening proxy for the
                       # strength-only naive25 sizer: the sizer runs with
                       # deflection_limit=false, and apply_naive_drf applies
                       # δ/NAIVE_DRF vs L/360 (and L/240) here.

# ── Run sizing pipelines ─────────────────────────────────────────────────────

mkpath(results_path)

new_rows = DataFrame(
    name            = String[],
    variant         = String[],
    max_depth       = Float64[],
    collinear       = Bool[],
    beam_sizer      = String[],
    steel_norm      = Float64[],
    concrete_norm   = Float64[],
    rebar_norm      = Float64[],
    column_norm     = Float64[],
    area            = Float64[],
    n_beams         = Int[],
    composite       = Bool[],
    staged_conv     = Bool[],
    n_L360_fail     = Int[],
    n_L240_fail     = Int[],
    max_util_M      = Float64[],
    max_util_V      = Float64[],
    max_δ_total_mm  = Float64[],
    global_δ_ok     = Bool[],
    result_ok       = Bool[],
)

"""
    apply_naive_drf(r, drf)

Return a NamedTuple of pass/fail metrics for a post-processed result `r`,
assuming the computed deflections are reduced by `drf` as a naive composite-
stiffening proxy. The underlying sized sections are unchanged — only the
serviceability verdict is recomputed.

The deflection check is applied unconditionally here: the naive path sizes
with `deflection_limit=false`, so `r.deflection_limit` is `false` and the
normal postprocess serviceability gate is bypassed by design. This helper
reinstates that check with the δ/`drf` reduction.
"""
function apply_naive_drf(r, drf::Real)
    δ_live_adj  = r.δ_live  ./ drf
    δ_total_adj = r.δ_total ./ drf

    live_ok  = δ_live_adj  .<= r.Δ_limit_live
    total_ok = δ_total_adj .<= r.Δ_limit_total
    n_L360_fail = count(.!live_ok)
    n_L240_fail = count(.!total_ok)

    max_δ_total = isempty(δ_total_adj) ? 0.0 : maximum(δ_total_adj)
    global_δ_ok = r.max_bay_span <= 0 || max_δ_total <= r.max_bay_span / 180.0

    strength_ok       = r.max_util_M <= 1.0 + 1e-3 && r.max_util_V <= 1.0 + 1e-3
    column_ok         = r.max_col_util <= 1.0 + 1e-3
    serviceability_ok = n_L360_fail == 0 && n_L240_fail == 0 && global_δ_ok
    result_ok = r.area > 0 && !isempty(r.ids) &&
                strength_ok && column_ok && serviceability_ok

    return (; n_L360_fail, n_L240_fail, max_δ_total, global_δ_ok, result_ok)
end

total_configs = length(geom_files) * length(max_depths) * length(variants)
global done_count = 0

for json_file in geom_files
    path = joinpath(main_path, json_file)
    name = replace(json_file, ".json" => "")

    json_string = replace(read(path, String), "\\n" => "")
    geometry_dict = JSON.parse(JSON.parse(json_string))
    if !(geometry_dict isa Dict)
        geometry_dict = convert(Dict{String,Any}, geometry_dict)
    end

    for max_depth in max_depths, variant in variants
        global done_count += 1
        println("\n[$done_count/$total_configs] $name — max_depth=$(max_depth)in — $(variant.name)")

        try
            geom, _type_info = generate_from_json(geometry_dict; plot=false, drawn=false)

            slab_params = SlabAnalysisParams(
                geom,
                slab_name       = name,
                slab_type       = slab_type,
                vector_1d       = vector_1d,
                slab_sizer      = slab_sizer,
                spacing         = 0.1,
                plot_analysis   = false,
                fix_param       = true,
                slab_units      = :m,
            )

            sizing_params = SlabSizingParams(
                live_load                   = psf_to_ksi(50),
                superimposed_dead_load      = psf_to_ksi(15),
                slab_dead_load              = 0.0,
                live_factor                 = 1.6,
                dead_factor                 = 1.2,
                beam_sizer                  = :continuous,
                nlp_solver                  = :Ipopt,
                max_depth                   = max_depth,
                beam_units                  = :in,
                serviceability_lim          = 360,
                collinear                   = false,
                minimum_continuous          = true,
                n_max_sections              = 0,
                composite_action            = variant.composite,
                deflection_reduction_factor = variant.drf,
                deflection_limit            = variant.deflection_limit,
            )

            t0 = time()
            slab_params = analyze_slab(slab_params)
            slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)
            dt_size = time() - t0

            if isempty(sizing_params.minimizers)
                println("  ⚠ No feasible sections found — skipping")
                continue
            end

            res_noncol = postprocess_slab(slab_params, sizing_params, check_collinear=false)
            res_col    = postprocess_slab(slab_params, sizing_params, check_collinear=true)

            # Persist raw (strict post-processing) results in both variants — the
            # naive25 override is applied in-memory only for the comparison table.
            append_results_to_csv(results_path,
                results_name * "_" * variant.name,
                [res_noncol, res_col])

            for (r, coll) in [(res_noncol, false), (res_col, true)]
                if variant.name == "naive25"
                    adj = apply_naive_drf(r, NAIVE_DRF)
                    max_δ_total_mm = adj.max_δ_total * 25.4
                    n_L360_fail    = adj.n_L360_fail
                    n_L240_fail    = adj.n_L240_fail
                    global_δ_ok    = adj.global_δ_ok
                    result_ok      = adj.result_ok
                else
                    max_δ_total_mm = isempty(r.δ_total) ? 0.0 : maximum(r.δ_total) * 25.4
                    n_L360_fail    = r.n_L360_fail
                    n_L240_fail    = r.n_L240_fail
                    global_δ_ok    = r.global_δ_ok
                    result_ok      = r.result_ok
                end

                push!(new_rows, (
                    name, variant.name, Float64(max_depth), coll, String(r.beam_sizer),
                    r.norm_mass_beams, r.norm_mass_slab, r.norm_mass_rebar,
                    r.norm_mass_columns, r.area, length(r.ids),
                    r.composite_action, r.staged_converged,
                    n_L360_fail, n_L240_fail,
                    r.max_util_M, r.max_util_V,
                    round(max_δ_total_mm, digits=2),
                    global_δ_ok, result_ok,
                ))
            end

            println("  ✓ sized in $(round(dt_size, digits=1))s — " *
                    "steel=$(round(res_col.norm_mass_beams, digits=2)) kg/m² (collinear, $(variant.name))")

        catch e
            @warn "Failed" name=name max_depth=max_depth variant=variant.name exception=(e, catch_backtrace())
        end

        GC.gc()
    end
end

# ── Load old baseline ────────────────────────────────────────────────────────

println("\n\n══════════════════════════════════════════════")
println("  Loading legacy dev-remote baseline (bare steel, strict L/360)")
println("══════════════════════════════════════════════\n")

old_df = CSV.read(old_csv_path, DataFrame)

old_filtered = filter(row ->
    row.slab_type == "isotropic" &&
    row.slab_sizer == "uniform" &&
    row.beam_sizer == "continuous" &&
    row.vector_1d_x == 0.0 &&
    row.vector_1d_y == 0.0 &&
    row.collinear == true &&
    row.area > 0 &&
    row.steel_norm > 0,
    old_df
)

target_names = Set(replace(f, ".json" => "") for f in geom_files)
old_subset = filter(row -> row.name in target_names, old_filtered)

println("Old baseline: $(nrow(old_subset)) rows for $(length(target_names)) geometries")

# ── Build comparison table ───────────────────────────────────────────────────

new_col = filter(row -> row.collinear == true, new_rows)

comp_rows   = filter(row -> row.variant == "composite", new_col)
naive_rows  = filter(row -> row.variant == "naive25",   new_col)

comparison = DataFrame(
    geometry         = String[],
    max_depth        = Int[],
    old_steel        = Float64[],
    naive_steel      = Float64[],
    comp_steel       = Float64[],
    Δ_comp_vs_old_pct = Float64[],
    Δ_comp_vs_naive_pct = Float64[],
    old_concrete     = Float64[],
    comp_concrete    = Float64[],
    comp_n_beams     = Int[],
    staged_conv      = Bool[],
    comp_L360_fail   = Int[],
    comp_L240_fail   = Int[],
    comp_max_δ_mm    = Float64[],
    naive_L360_fail  = Int[],
    naive_L240_fail  = Int[],
    naive_max_δ_mm   = Float64[],
    comp_ok          = Bool[],
    naive_ok         = Bool[],
)

for comp_r in eachrow(comp_rows)
    naive_match = filter(row ->
        row.name == comp_r.name && row.max_depth == comp_r.max_depth,
        naive_rows,
    )
    old_match = filter(row ->
        row.name == comp_r.name && row.max_depth == comp_r.max_depth,
        old_subset,
    )

    old_steel = nrow(old_match) == 0 ? NaN : old_match[1, :steel_norm]
    old_conc  = (nrow(old_match) == 0 || !hasproperty(old_match, :concrete_norm)) ?
                NaN : old_match[1, :concrete_norm]

    naive_steel     = nrow(naive_match) == 0 ? NaN : naive_match[1, :steel_norm]
    naive_L360      = nrow(naive_match) == 0 ? -1  : naive_match[1, :n_L360_fail]
    naive_L240      = nrow(naive_match) == 0 ? -1  : naive_match[1, :n_L240_fail]
    naive_max_δ     = nrow(naive_match) == 0 ? NaN : naive_match[1, :max_δ_total_mm]
    naive_ok        = nrow(naive_match) == 0 ? false : naive_match[1, :result_ok]

    Δ_comp_old = (isnan(old_steel) || old_steel == 0) ? NaN :
        round((comp_r.steel_norm - old_steel) / old_steel * 100, digits=1)
    Δ_comp_naive = (isnan(naive_steel) || naive_steel == 0) ? NaN :
        round((comp_r.steel_norm - naive_steel) / naive_steel * 100, digits=1)

    push!(comparison, (
        comp_r.name, Int(comp_r.max_depth),
        round(old_steel, digits=2),
        round(naive_steel, digits=2),
        round(comp_r.steel_norm, digits=2),
        Δ_comp_old, Δ_comp_naive,
        round(old_conc, digits=1),
        round(comp_r.concrete_norm, digits=1),
        comp_r.n_beams, comp_r.staged_conv,
        comp_r.n_L360_fail, comp_r.n_L240_fail, comp_r.max_δ_total_mm,
        naive_L360, naive_L240, naive_max_δ,
        comp_r.result_ok, naive_ok,
    ))
end

# ── Print results ────────────────────────────────────────────────────────────

println("\n\n╔══════════════════════════════════════════════════════════════════╗")
println("║  COMPOSITE vs NAIVE-2.5× (new pipeline) vs legacy L/360 (old CSV)║")
println("║                  Collinear • Continuous                          ║")
println("╚══════════════════════════════════════════════════════════════════╝\n")

println("Composite:  composite_action=true,  drf=1.0 (staged L/360 + L/240)")
println("Naive-2.5×: bare steel, deflection_limit=OFF in sizing;")
println("            post-processing: δ/2.5 vs L/360 & L/240")
println("Old CSV:    bare steel, strict L/360 on bare section (dev-remote legacy)\n")

if nrow(comparison) > 0
    pretty_table(comparison;
        column_labels = ["Geom", "Depth",
                         "Old steel", "Naive steel", "Comp steel",
                         "Δ comp vs old%", "Δ comp vs naive%",
                         "Old conc", "Comp conc", "#Beams", "Staged",
                         "Comp L360", "Comp L240", "Comp δ mm",
                         "Naive L360", "Naive L240", "Naive δ mm",
                         "Comp OK", "Naive OK"],
        alignment = [:l, :r, :r, :r, :r, :r, :r, :r, :r, :r, :c,
                     :r, :r, :r, :r, :r, :r, :c, :c],
    )

    valid_old = filter(row -> !isnan(row.Δ_comp_vs_old_pct),   comparison)
    valid_nv  = filter(row -> !isnan(row.Δ_comp_vs_naive_pct), comparison)
    if nrow(valid_old) > 0
        println("\nComposite vs legacy bare-L/360 CSV ($(nrow(valid_old)) configs):")
        println("  Mean  Δ steel: $(round(mean(valid_old.Δ_comp_vs_old_pct),   digits=1))%")
        println("  Median Δ steel: $(round(median(valid_old.Δ_comp_vs_old_pct), digits=1))%")
    end
    if nrow(valid_nv) > 0
        println("\nComposite vs naive-2.5× (new) ($(nrow(valid_nv)) configs):")
        println("  Mean  Δ steel: $(round(mean(valid_nv.Δ_comp_vs_naive_pct),   digits=1))%")
        println("  Median Δ steel: $(round(median(valid_nv.Δ_comp_vs_naive_pct), digits=1))%")
    end

    println("\nServiceability (composite):")
    println("  L/360 fail: $(count(comparison.comp_L360_fail .> 0))/$(nrow(comparison))")
    println("  L/240 fail: $(count(comparison.comp_L240_fail .> 0))/$(nrow(comparison))")
    println("  result_ok:  $(count(comparison.comp_ok))/$(nrow(comparison))")

    println("\nServiceability (naive-2.5×, δ/2.5 pass/fail):")
    naive_valid = filter(row -> row.naive_L360_fail >= 0, comparison)
    if nrow(naive_valid) > 0
        println("  L/360 fail: $(count(naive_valid.naive_L360_fail .> 0))/$(nrow(naive_valid))")
        println("  L/240 fail: $(count(naive_valid.naive_L240_fail .> 0))/$(nrow(naive_valid))")
        println("  result_ok:  $(count(naive_valid.naive_ok))/$(nrow(naive_valid))")
    else
        println("  (no naive-2.5× runs succeeded)")
    end
else
    println("No comparison rows generated — check that geometries ran successfully.")
end

# ── Save comparison CSV ──────────────────────────────────────────────────────
comparison_path = joinpath(results_path, "comparison_summary.csv")
CSV.write(comparison_path, comparison)
new_rows_path = joinpath(results_path, "new_rows_all_variants.csv")
CSV.write(new_rows_path, new_rows)
println("\nComparison saved to: $comparison_path")
println("Per-variant/per-collinearity rows saved to: $new_rows_path")
println("Full composite results saved to: $(joinpath(results_path, results_name * "_composite.csv"))")
println("Full naive-2.5× results saved to: $(joinpath(results_path, results_name * "_naive25.csv"))")
println("\nDone.")
