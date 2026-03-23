"""
Run topology-aware beam optimization (`optimize_indeterminate`) for one or more
JSON slab layouts and save concise outputs.

Outputs (per geometry):
- `*_topopt_beams.csv`: optimized section geometry and per-beam demand summary

Global output:
- `topopt_summary.csv`: one-row summary per geometry
"""

include("_scripts.jl")

# ── configuration ─────────────────────────────────────────────────────────────
main_path = "Geometries/topology/"
results_path = "SlabDesignFactors/results/topopt_runs/"
mkpath(results_path)

run_all_geometries = false
single_geometry = "r1c2.json"
# If provided, this takes precedence over `main_path` + `single_geometry`.
# Example:
# single_geometry_path = "Geometries/special/topopt_test2.json"
single_geometry_path = "Geometries/special/topopt_test2.json"

slab_type = :isotropic
vector_1d = [1.0, 0.0]
slab_sizer = :uniform
spacing = 0.1

max_depth_in = 40.0
serviceability_lim = 360.0

# Loads (ksi)
live_load = psf_to_ksi(50)
superimposed_dead_load = psf_to_ksi(15)
slab_dead_load = 0.0
live_factor = 1.6
dead_factor = 1.2

composite_action = false
deflection_reduction_factor = 1.0

# Continuation schedule for topopt. Start from default and tweak as needed.
continuation_schedule = topopt_default_schedule()
# continuation_schedule[3] = merge(continuation_schedule[3], (penalty=8.0, maxeval=3000))

# ── helpers ───────────────────────────────────────────────────────────────────
function load_geometry_from_json(path::String)
    raw = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))
    geometry_dict = raw isa Dict ? raw : Dict(pairs(raw))
    geometry, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)
    return geometry
end

function build_slab_params(geometry, name::String)
    return SlabAnalysisParams(
        geometry,
        slab_name=name,
        slab_type=slab_type,
        vector_1d=vector_1d,
        slab_sizer=slab_sizer,
        spacing=spacing,
        plot_analysis=false,
        fix_param=true,
        slab_units=:m,
    )
end

function build_sizing_params()
    return SlabSizingParams(
        live_load=live_load,
        superimposed_dead_load=superimposed_dead_load,
        slab_dead_load=slab_dead_load,
        live_factor=live_factor,
        dead_factor=dead_factor,
        beam_sizer=:continuous, # topopt routine is independent from discrete/continuous toggle
        max_depth=max_depth_in,
        beam_units=:in,
        serviceability_lim=serviceability_lim,
        collinear=false,
        minimum_continuous=true,
        n_max_sections=0,
        composite_action=composite_action,
        deflection_reduction_factor=deflection_reduction_factor,
    )
end

# ── run list ──────────────────────────────────────────────────────────────────
runs = if run_all_geometries
    jsons = filter(x -> endswith(x, ".json"), readdir(main_path))
    [(name=replace(f, ".json" => ""), path=joinpath(main_path, f)) for f in jsons]
elseif !isempty(single_geometry_path)
    [(name=replace(basename(single_geometry_path), ".json" => ""), path=single_geometry_path)]
else
    [(name=replace(single_geometry, ".json" => ""), path=joinpath(main_path, single_geometry))]
end

if isempty(runs)
    error("No JSON files found to run (main_path=$main_path)")
end

summary_df = DataFrame(
    geometry=String[],
    elapsed_s=Float64[],
    n_beams=Int[],
    total_volume_in3=Float64[],
    mean_area_in2=Float64[],
    max_area_in2=Float64[],
    min_area_in2=Float64[],
    max_Mu=Float64[],
    max_Vu=Float64[],
)

for run_cfg in runs
    name = run_cfg.name
    path = run_cfg.path

    println("\n================================================")
    println("TOPOPT RUN: $name")
    println("================================================")

    if !isfile(path)
        @warn "Skipping missing geometry file: $path"
        continue
    end

    t0 = time()
    try
        geometry = load_geometry_from_json(path)
        slab_params = build_slab_params(geometry, name)
        sizing_params = build_sizing_params()

        minimizers, sizing_params = optimize_indeterminate(
            slab_params,
            sizing_params;
            continuation_schedule=continuation_schedule,
        )

        elapsed_s = time() - t0
        areas = [I_symm(m...).A for m in minimizers]
        total_volume_in3 = sum(sizing_params.minimums)
        n_beams = length(minimizers)

        beam_df = DataFrame(
            beam_idx=1:n_beams,
            h=[m[1] for m in minimizers],
            w=[m[2] for m in minimizers],
            tw=[m[3] for m in minimizers],
            tf=[m[4] for m in minimizers],
            A_in2=areas,
            volume_in3=sizing_params.minimums,
            section_id=sizing_params.ids,
            M_max=sizing_params.M_maxs,
            V_max=sizing_params.V_maxs,
            x_max=sizing_params.x_maxs,
        )

        beam_csv = joinpath(results_path, "$(name)_topopt_beams.csv")
        CSV.write(beam_csv, beam_df)

        push!(summary_df, (
            name,
            elapsed_s,
            n_beams,
            total_volume_in3,
            mean(areas),
            maximum(areas),
            minimum(areas),
            isempty(sizing_params.M_maxs) ? NaN : maximum(sizing_params.M_maxs),
            isempty(sizing_params.V_maxs) ? NaN : maximum(sizing_params.V_maxs),
        ))

        println("✓ Completed in $(round(elapsed_s, digits=2)) s")
        println("  n_beams=$n_beams, total_volume_in3=$(round(total_volume_in3, digits=2))")
        println("  Saved: $beam_csv")

    catch e
        elapsed_s = time() - t0
        @warn "Topopt run failed for $name after $(round(elapsed_s, digits=2)) s" exception=(e, catch_backtrace())
    end
end

summary_csv = joinpath(results_path, "topopt_summary.csv")
CSV.write(summary_csv, summary_df)
println("\n✓ Topopt summary saved to $summary_csv")
