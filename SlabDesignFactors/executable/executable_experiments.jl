"""
Cluster-ready experiment runner.

Experiments included:
1. `max_depths` (existing study retained)
2. `strip_resolution` sensitivity
3. `constrained_inventory` sensitivities:
   - `constrained_inventory_resized` (n_max_sections sweep)
   - `constrained_inventory_max_depth` (max_depth sweep)
4. `nlp_solver_comparison` (continuous sizing across NLP solvers with MIP warm-start)
5. `material_scenario_mc` (concrete material scenarios × slab sizers with
   per-material Monte Carlo on ECCs for steel, concrete, and rebar)
6. `validation_mip` (validation JSONs with drawn=true, perimeter beams, holes,
   element IDs, building-type loads, and unconstrained MIP sizing)

All studies are resumable from CSV outputs and safe for repeated 12-hour
Slurm re-submissions.
"""
include("../SlabDesignFactors.jl")

using Base.Threads
using LinearAlgebra
using JSON
using CSV
using DataFrames
using Statistics
using Random
using .SlabDesignFactors

const DEFAULT_TOPOLOGY_JSON_PATH = "/home/nhirt/2024_Slab-Design-Factors/SlabDesignFactors/jsons/topology/"
const DEFAULT_VALIDATION_JSON_PATH = "/home/nhirt/2024_Slab-Design-Factors/Geometries/validation/"
const DEFAULT_STUDIES = Set(["max_depths", "strip_resolution", "constrained_inventory", "nlp_solver_comparison", "material_scenario_mc", "validation_mip"])
const NLP_SOLVER_LIST = [:MIP, :MMA, :SLSQP, :CCSAQ, :COBYLA, :Ipopt]

# ── Shared experiment defaults (topology sweeps; see `full_sweep_defaults.jl`) ─
const DEFAULT_LIVE_LOAD   = SlabDesignFactors.FULL_SWEEP_LIVE_LOAD
const DEFAULT_SDL         = SlabDesignFactors.FULL_SWEEP_SUPERIMPOSED_DEAD_LOAD
const DEFAULT_LIVE_FACTOR = SlabDesignFactors.FULL_SWEEP_LIVE_FACTOR
const DEFAULT_DEAD_FACTOR = SlabDesignFactors.FULL_SWEEP_DEAD_FACTOR
const DEFAULT_MAX_DEPTH   = 40.0
const DEFAULT_SERV_LIM    = SlabDesignFactors.FULL_SWEEP_SERVICEABILITY_LIM
const DEFAULT_SLAB_TYPE   = :isotropic
const DEFAULT_VECTOR_1D   = [0.0, 0.0]
const DEFAULT_SPACING     = SlabDesignFactors.FULL_SWEEP_STRIP_SPACING
const DEFAULT_COMPOSITE   = SlabDesignFactors.FULL_SWEEP_COMPOSITE_ACTION
const DEFAULT_COLLINEAR   = SlabDesignFactors.FULL_SWEEP_COLLINEAR
const DEFAULT_DRF         = SlabDesignFactors.FULL_SWEEP_DEFLECTION_REDUCTION_FACTOR

# Columns added to custom study CSVs (strip / NLP / material) for parity with `create_results_dataframe`.
const STAGING_SUMMARY_COLS = [
    :beam_sizer, :nlp_solver, :deflection_limit, :composite_action, :staged_converged,
    :staged_n_violations, :n_L360_fail, :n_L240_fail, :global_δ_ok, :max_staged_δ_total_mm,
]

# ── Monte Carlo parameters ──────────────────────────────────────────────────
const MC_N_SAMPLES   = 10_000
const MC_CV_STEEL    = 0.10
const MC_CV_REBAR    = 0.10
const MC_CV_CONCRETE = 0.15
const MC_ECC_STEEL_NOM = SlabDesignFactors.ECC_STEEL  # 1.22  kgCO₂e/kg
const MC_ECC_REBAR_NOM = SlabDesignFactors.ECC_REBAR  # 0.854 kgCO₂e/kg

"""
    ensure_dir(path::String)

Create output directory if missing.
"""
function ensure_dir(path::String)
    if !isdir(path)
        mkpath(path)
    end
end

"""
    read_geometry_from_json(path::String) -> Dict

Read and parse the geometry JSON payload format used in this repository.
"""
function read_geometry_from_json(path::String)
    json_string = replace(read(path, String), "\\n" => "")
    return JSON.parse(JSON.parse(json_string, dicttype=Dict))
end

"""
    parse_studies_arg(raw::String) -> Set{String}

Parse an optional comma-separated studies argument.
"""
function parse_studies_arg(raw::String)::Set{String}
    studies = Set{String}()
    for item in split(raw, ",")
        s = strip(lowercase(item))
        isempty(s) && continue
        push!(studies, s)
    end
    return studies
end

"""
    write_table_atomically(path, df)

Write a DataFrame to CSV via temporary file then rename.
"""
function write_table_atomically(path::String, df::DataFrame)
    tmp = path * ".tmp"
    CSV.write(tmp, df)
    Base.Filesystem.rename(tmp, path)
end

"""
    ensure_results_columns!(df, colnames)

Append missing columns to `df` when upgrading CSV schemas so existing cluster
outputs remain readable. Existing rows get `missing` in new columns.
"""
function ensure_results_columns!(df::DataFrame, colnames::Vector{Symbol})
    for sym in colnames
        if !hasproperty(df, sym)
            df[!, sym] = [missing for _ in 1:nrow(df)]
        end
    end
    nothing
end

function done_set_max_depth(results_file::String)::Set{Tuple{String, String, Float64, String, Float64, Float64}}
    out = Set{Tuple{String, String, Float64, String, Float64, Float64}}()
    if !isfile(results_file)
        return out
    end
    try
        df = CSV.read(results_file, DataFrame)
        for r in eachrow(df)
            push!(out, (
                String(r.name),
                String(r.slab_sizer),
                Float64(r.max_depth),
                String(r.slab_type),
                Float64(r.vector_1d_x),
                Float64(r.vector_1d_y),
            ))
        end
    catch e
        @warn "Failed to read max_depth done-set; rerunning from scratch for this study" file=results_file exception=e
    end
    return out
end

function done_set_constrained(results_file::String)::Set{Tuple{String, Int, Int}}
    out = Set{Tuple{String, Int, Int}}()
    if !isfile(results_file)
        return out
    end
    try
        df = CSV.read(results_file, DataFrame)
        if !all(c -> c in names(df), [:name, :max_depth, :unique_sections])
            return out
        end
        for r in eachrow(df)
            push!(out, (String(r.name), Int(round(r.max_depth)), Int(round(r.unique_sections))))
        end
    catch e
        @warn "Failed to read constrained inventory done-set; rerunning from scratch for this study" file=results_file exception=e
    end
    return out
end

function done_set_strip_resolution(results_file::String)::Set{Tuple{String, Float64}}
    out = Set{Tuple{String, Float64}}()
    if !isfile(results_file)
        return out
    end
    try
        df = CSV.read(results_file, DataFrame)
        if !all(c -> c in names(df), [:geometry, :spacing])
            return out
        end
        for r in eachrow(df)
            push!(out, (String(r.geometry), Float64(r.spacing)))
        end
    catch e
        @warn "Failed to read strip-resolution done-set; rerunning from scratch for this study" file=results_file exception=e
    end
    return out
end

function done_set_nlp_solver(results_file::String)::Set{Tuple{String, String, String, String, Float64}}
    out = Set{Tuple{String, String, String, String, Float64}}()
    if !isfile(results_file)
        return out
    end
    try
        df = CSV.read(results_file, DataFrame)
        required = [:geometry, :solver, :slab_type, :slab_sizer, :max_depth]
        if !all(c -> c in names(df), required)
            return out
        end
        for r in eachrow(df)
            push!(out, (
                String(r.geometry),
                String(r.solver),
                String(r.slab_type),
                String(r.slab_sizer),
                Float64(r.max_depth),
            ))
        end
    catch e
        @warn "Failed to read NLP solver comparison done-set; rerunning from scratch for this study" file=results_file exception=e
    end
    return out
end

function run_max_depths(results_path::String; json_path::String=DEFAULT_TOPOLOGY_JSON_PATH)
    println("\n=== Study: max_depths ===")
    ensure_dir(results_path)

    results_name = "max_depths"
    results_file = joinpath(results_path, results_name * ".csv")
    done = done_set_max_depth(results_file)

    slab_types = [:isotropic]
    vector_1ds = [[0.0, 0.0]]
    max_depths = [10, 15, 20, 25, 30, 35, 40, 45, 50, 10000]
    slab_sizers = [:uniform, :cellular]

    configs = NamedTuple[]
    sub_paths = sort(filter(x -> endswith(x, ".json"), readdir(json_path)))
    for sub_path in sub_paths
        path = joinpath(json_path, sub_path)
        name = replace(sub_path, ".json" => "")
        for max_depth in max_depths
            for slab_sizer in slab_sizers
                for (i, slab_type) in enumerate(slab_types)
                    vector_1d = vector_1ds[i]
                    key = (name, String(slab_sizer), Float64(max_depth), String(slab_type), Float64(vector_1d[1]), Float64(vector_1d[2]))
                    if key in done
                        continue
                    end
                    push!(configs, (
                        path=path,
                        name=name,
                        slab_type=slab_type,
                        vector_1d=vector_1d,
                        slab_sizer=slab_sizer,
                        max_depth=max_depth,
                    ))
                end
            end
        end
    end

    println("Pending max_depth configs: $(length(configs))")
    isempty(configs) && return

    append_lock = ReentrantLock()
    Threads.@threads for i in eachindex(configs)
        cfg = configs[i]
        println("[max_depths] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name) $(cfg.slab_type) $(cfg.slab_sizer) $(cfg.max_depth)in")
        try
            geometry_dict = read_geometry_from_json(cfg.path)
            geometry = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=false)

            slab_params = SlabDesignFactors.SlabAnalysisParams(
                geometry,
                slab_name=cfg.name,
                slab_type=cfg.slab_type,
                vector_1d=cfg.vector_1d,
                slab_sizer=cfg.slab_sizer,
                spacing=DEFAULT_SPACING,
                plot_analysis=false,
                fix_param=true,
                slab_units=:m,
            )

            beam_sizing_params = SlabDesignFactors.SlabSizingParams(
                live_load=DEFAULT_LIVE_LOAD,
                superimposed_dead_load=DEFAULT_SDL,
                slab_dead_load=0.0,
                live_factor=DEFAULT_LIVE_FACTOR,
                dead_factor=DEFAULT_DEAD_FACTOR,
                beam_sizer=:discrete,
                max_depth=cfg.max_depth,
                beam_units=:in,
                serviceability_lim=DEFAULT_SERV_LIM,
                minimum_continuous=true,
                collinear=DEFAULT_COLLINEAR,
                composite_action=DEFAULT_COMPOSITE,
                deflection_reduction_factor=DEFAULT_DRF,
            )

            iteration_result = collect(SlabDesignFactors.iterate_discrete_continuous(slab_params, beam_sizing_params))
            lock(append_lock) do
                SlabDesignFactors.append_results_to_csv(results_path * "/", results_name, iteration_result)
            end
        catch e
            @warn "Failed max_depth config" name=cfg.name slab_type=cfg.slab_type slab_sizer=cfg.slab_sizer max_depth=cfg.max_depth exception=e
        end
    end
end

function run_strip_resolution(results_path::String; json_path::String=DEFAULT_TOPOLOGY_JSON_PATH)
    println("\n=== Study: strip_resolution ===")
    ensure_dir(results_path)

    results_file = joinpath(results_path, "strip_resolution_sensitivity.csv")
    done = done_set_strip_resolution(results_file)
    spacings = [1.0, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001]

    # Seed existing table for lock-protected, atomic rewrites.
    out_df = if isfile(results_file)
        CSV.read(results_file, DataFrame)
    else
        DataFrame(
            geometry=String[],
            spacing=Float64[],
            n_cells=Int[],
            n_beams=Int[],
            n_strips=Int[],
            mass_beams_kg=Float64[],
            norm_mass_beams=Float64[],
            ec_steel=Float64[],
            ec_slab=Float64[],
            ec_rebar=Float64[],
            ec_total=Float64[],
            max_deflection_in=Float64[],
            slab_area_m2=Float64[],
            mean_slab_depth_m=Float64[],
            elapsed_s=Float64[],
        )
    end
    ensure_results_columns!(out_df, STAGING_SUMMARY_COLS)

    configs = NamedTuple[]
    sub_paths = sort(filter(x -> endswith(x, ".json"), readdir(json_path)))
    for sub_path in sub_paths
        name = replace(sub_path, ".json" => "")
        path = joinpath(json_path, sub_path)
        for spacing in spacings
            key = (name, spacing)
            if key in done
                continue
            end
            push!(configs, (name=name, path=path, spacing=spacing))
        end
    end

    println("Pending strip-resolution configs: $(length(configs))")
    isempty(configs) && return

    table_lock = ReentrantLock()
    Threads.@threads for i in eachindex(configs)
        cfg = configs[i]
        t0 = time()
        println("[strip_resolution] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name) spacing=$(cfg.spacing)")

        try
            geometry_dict = read_geometry_from_json(cfg.path)
            geometry = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=false)
            slab_params = SlabDesignFactors.SlabAnalysisParams(
                geometry,
                slab_name=cfg.name,
                slab_type=:isotropic,
                vector_1d=[1.0, 0.0],
                slab_sizer=:uniform,
                spacing=cfg.spacing,
                plot_analysis=false,
                fix_param=true,
                slab_units=:m,
            )

            sizing_params = SlabDesignFactors.SlabSizingParams(
                live_load=DEFAULT_LIVE_LOAD,
                superimposed_dead_load=DEFAULT_SDL,
                slab_dead_load=0.0,
                live_factor=DEFAULT_LIVE_FACTOR,
                dead_factor=DEFAULT_DEAD_FACTOR,
                beam_sizer=:discrete,
                max_depth=DEFAULT_MAX_DEPTH,
                beam_units=:in,
                serviceability_lim=DEFAULT_SERV_LIM,
                collinear=true,
                minimum_continuous=true,
                n_max_sections=0,
                composite_action=true,
                deflection_reduction_factor=1.0,
            )

            slab_params = SlabDesignFactors.analyze_slab(slab_params)
            slab_params, sizing_params = SlabDesignFactors.optimal_beamsizer(slab_params, sizing_params)
            results = SlabDesignFactors.postprocess_slab(slab_params, sizing_params, check_collinear=false)

            n_beams = length(results.Δ_local)
            max_δ = maximum(maximum(abs.(d)) for d in results.Δ_local if !isempty(d); init=0.0)
            n_strips = length(slab_params.load_areas)
            mean_depth = isempty(slab_params.slab_depths) ? 0.0 : mean(slab_params.slab_depths)
            dt = time() - t0

            max_staged_mm = round(results.max_δ_total * 25.4, digits=3)
            row = (
                cfg.name,
                cfg.spacing,
                0, # kept for compatibility with plotting scripts; this runner does not precompute cells.
                n_beams,
                n_strips,
                results.mass_beams,
                results.norm_mass_beams,
                results.embodied_carbon_beams,
                results.embodied_carbon_slab,
                results.embodied_carbon_rebar,
                results.embodied_carbon_beams + results.embodied_carbon_slab + results.embodied_carbon_rebar,
                max_δ,
                results.area,
                mean_depth,
                dt,
                String(results.beam_sizer),
                results.nlp_solver,
                results.deflection_limit,
                results.composite_action,
                results.staged_converged,
                results.staged_n_violations,
                results.n_L360_fail,
                results.n_L240_fail,
                results.global_δ_ok,
                max_staged_mm,
            )

            lock(table_lock) do
                push!(out_df, row)
                write_table_atomically(results_file, out_df)
            end
        catch e
            dt = time() - t0
            lock(table_lock) do
                push!(out_df, (
                    cfg.name, cfg.spacing, 0, 0, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, dt,
                    missing, missing, missing, missing, missing, missing, missing, missing, missing, missing,
                ))
                write_table_atomically(results_file, out_df)
            end
            @warn "Failed strip-resolution config" name=cfg.name spacing=cfg.spacing exception=e
        end
    end
end

function run_constrained_inventory(results_path::String; json_path::String=DEFAULT_VALIDATION_JSON_PATH)
    println("\n=== Study: constrained_inventory ===")
    ensure_dir(results_path)

    live_load_dict = Dict(
        "office" => 60,
        "school" => 50,
        "warehouse" => 250,
    )
    names = ["office", "school", "warehouse"]

    # Substudy A: n_max_sections sweep.
    resized_results_name = "constrained_inventory_resized"
    resized_file = joinpath(results_path, resized_results_name * ".csv")
    done_resized = done_set_constrained(resized_file)
    resized_configs = NamedTuple[]
    for name in names
        path = joinpath(json_path, name * ".json")
        for n_sections in 0:30
            key = (name, 40, n_sections)
            if !(key in done_resized)
                push!(resized_configs, (name=name, path=path, max_depth=40, n_max_sections=n_sections))
            end
        end
    end

    # Substudy B: max_depth sweep.
    max_depth_results_name = "constrained_inventory_max_depth"
    max_depth_file = joinpath(results_path, max_depth_results_name * ".csv")
    done_max_depth = done_set_constrained(max_depth_file)
    max_depth_values = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 1000]
    max_depth_configs = NamedTuple[]
    for name in names
        path = joinpath(json_path, name * ".json")
        for depth in max_depth_values
            key = (name, depth, depth)
            if !(key in done_max_depth)
                push!(max_depth_configs, (name=name, path=path, max_depth=depth, n_max_sections=0))
            end
        end
    end

    println("Pending constrained_inventory resized configs: $(length(resized_configs))")
    println("Pending constrained_inventory max_depth configs: $(length(max_depth_configs))")

    append_lock = ReentrantLock()

    # Helper closure shared by both substudies.
    function run_constrained_cfg(cfg::NamedTuple, results_name::String, unique_sections_value::Int)
        geometry_dict = read_geometry_from_json(cfg.path)
        geometry = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=false)

        slab_params = SlabDesignFactors.SlabAnalysisParams(
            geometry,
            slab_name=cfg.name,
            slab_type=:uniaxial,
            slab_thickness=4.75 * SlabDesignFactors.convert_to_m[:in],
            vector_1d=[0.0, 1.0],
            slab_sizer=:uniform,
            spacing=DEFAULT_SPACING,
            plot_analysis=false,
            fix_param=true,
            slab_units=:m,
        )

        beam_sizing_params = SlabDesignFactors.SlabSizingParams(
            live_load=SlabDesignFactors.psf_to_ksi(live_load_dict[cfg.name]),
            superimposed_dead_load=SlabDesignFactors.psf_to_ksi(20),
            slab_dead_load=SlabDesignFactors.psf_to_ksi(45),
            façade_load=SlabDesignFactors.plf_to_kpi(500),
            live_factor=1.6,
            dead_factor=1.2,
            beam_sizer=:discrete,
            max_depth=cfg.max_depth,
            beam_units=:in,
            serviceability_lim=360,
            collinear=false,
            minimum_continuous=true,
            n_max_sections=cfg.n_max_sections,
        )

        slab_params = SlabDesignFactors.analyze_slab(slab_params)
        slab_params, beam_sizing_params = SlabDesignFactors.optimal_beamsizer(slab_params, beam_sizing_params)
        if isempty(beam_sizing_params.minimizers)
            return
        end
        slab_results = SlabDesignFactors.postprocess_slab(slab_params, beam_sizing_params)
        lock(append_lock) do
            SlabDesignFactors.append_results_to_csv(results_path * "/", results_name, [slab_results], unique_sections=unique_sections_value)
        end
    end

    Threads.@threads for i in eachindex(resized_configs)
        cfg = resized_configs[i]
        println("[constrained_inventory:resized] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name) n_max_sections=$(cfg.n_max_sections)")
        try
            run_constrained_cfg(cfg, resized_results_name, cfg.n_max_sections)
        catch e
            @warn "Failed constrained inventory resized config" name=cfg.name n_max_sections=cfg.n_max_sections exception=e
        end
    end

    Threads.@threads for i in eachindex(max_depth_configs)
        cfg = max_depth_configs[i]
        println("[constrained_inventory:max_depth] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name) max_depth=$(cfg.max_depth)")
        try
            run_constrained_cfg(cfg, max_depth_results_name, cfg.max_depth)
        catch e
            @warn "Failed constrained inventory max-depth config" name=cfg.name max_depth=cfg.max_depth exception=e
        end
    end
end

function run_nlp_solver_comparison(results_path::String; json_path::String=DEFAULT_TOPOLOGY_JSON_PATH)
    println("\n=== Study: nlp_solver_comparison ===")
    ensure_dir(results_path)

    results_file = joinpath(results_path, "nlp_solver_comparison.csv")
    done = done_set_nlp_solver(results_file)

    # Keep physics/config fixed so solver differences are isolated.
    slab_type = :isotropic
    vector_1d = [1.0, 0.0]
    slab_sizer = :uniform
    spacing = DEFAULT_SPACING
    max_depth = 40.0

    out_df = if isfile(results_file)
        CSV.read(results_file, DataFrame)
    else
        DataFrame(
            geometry=String[],
            solver=String[],
            slab_type=String[],
            slab_sizer=String[],
            max_depth=Float64[],
            warm_start_used=Bool[],
            mip_seed_count=Int[],
            success=Bool[],
            elapsed_s=Float64[],
            n_beams=Int[],
            n_strips=Int[],
            norm_mass_beams=Float64[],
            ec_total=Float64[],
            max_deflection_in=Float64[],
            error_message=String[],
        )
    end
    ensure_results_columns!(out_df, STAGING_SUMMARY_COLS)

    configs = NamedTuple[]
    sub_paths = sort(filter(x -> endswith(x, ".json"), readdir(json_path)))
    for sub_path in sub_paths
        name = replace(sub_path, ".json" => "")
        path = joinpath(json_path, sub_path)
        for solver in NLP_SOLVER_LIST
            key = (name, String(solver), String(slab_type), String(slab_sizer), max_depth)
            if key in done
                continue
            end
            push!(configs, (name=name, path=path, solver=solver))
        end
    end

    println("Pending NLP solver comparison configs: $(length(configs))")
    isempty(configs) && return

    table_lock = ReentrantLock()
    Threads.@threads for i in eachindex(configs)
        cfg = configs[i]
        t0 = time()
        println("[nlp_solver_comparison] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name) solver=$(cfg.solver)")

        success = false
        warm_start_used = false
        mip_seed_count = 0
        n_beams = 0
        n_strips = 0
        norm_mass_beams = NaN
        ec_total = NaN
        max_deflection_in = NaN
        err = ""

        beam_sizer_str = missing
        nlp_solver_out = missing
        deflection_limit_out = missing
        composite_action_out = missing
        staged_converged_out = missing
        staged_n_violations_out = missing
        n_L360_fail_out = missing
        n_L240_fail_out = missing
        global_δ_ok_out = missing
        max_staged_mm_out = missing

        try
            geometry_dict = read_geometry_from_json(cfg.path)
            geometry = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=false)

            slab_params = SlabDesignFactors.SlabAnalysisParams(
                geometry,
                slab_name=cfg.name,
                slab_type=slab_type,
                vector_1d=vector_1d,
                slab_sizer=slab_sizer,
                spacing=spacing,
                plot_analysis=false,
                fix_param=true,
                slab_units=:m,
            )

            beam_sizer = cfg.solver == :MIP ? :discrete : :continuous
            nlp_solver = cfg.solver == :MIP ? :MMA : cfg.solver

            sizing_params = SlabDesignFactors.SlabSizingParams(
                live_load=DEFAULT_LIVE_LOAD,
                superimposed_dead_load=DEFAULT_SDL,
                slab_dead_load=0.0,
                live_factor=DEFAULT_LIVE_FACTOR,
                dead_factor=DEFAULT_DEAD_FACTOR,
                beam_sizer=beam_sizer,
                nlp_solver=nlp_solver,
                max_depth=max_depth,
                beam_units=:in,
                serviceability_lim=DEFAULT_SERV_LIM,
                minimum_continuous=true,
                collinear=true,
                composite_action=true,
                deflection_reduction_factor=1.0,
            )

            slab_params = SlabDesignFactors.analyze_slab(slab_params)
            slab_params, sizing_params = SlabDesignFactors.optimal_beamsizer(slab_params, sizing_params)

            if cfg.solver == :MIP
                warm_start_used = false
                mip_seed_count = 0
            else
                warm_start_used = !isnothing(sizing_params.mip_result) && !isempty(sizing_params.mip_result.ids)
                mip_seed_count = isnothing(sizing_params.mip_result) ? 0 : length(sizing_params.mip_result.ids)
            end

            if !isempty(sizing_params.minimizers)
                results = SlabDesignFactors.postprocess_slab(slab_params, sizing_params, check_collinear=false)
                n_beams = length(results.Δ_local)
                n_strips = length(slab_params.load_areas)
                norm_mass_beams = results.norm_mass_beams
                ec_total = results.embodied_carbon_beams + results.embodied_carbon_slab + results.embodied_carbon_rebar
                max_deflection_in = maximum(maximum(abs.(d)) for d in results.Δ_local if !isempty(d); init=0.0)
                success = true
                beam_sizer_str = String(results.beam_sizer)
                nlp_solver_out = results.nlp_solver
                deflection_limit_out = results.deflection_limit
                composite_action_out = results.composite_action
                staged_converged_out = results.staged_converged
                staged_n_violations_out = results.staged_n_violations
                n_L360_fail_out = results.n_L360_fail
                n_L240_fail_out = results.n_L240_fail
                global_δ_ok_out = results.global_δ_ok
                max_staged_mm_out = round(results.max_δ_total * 25.4, digits=3)
            else
                err = "No minimizers returned from optimizer"
            end
        catch e
            err = sprint(showerror, e)
            @warn "Failed NLP solver comparison config" name=cfg.name solver=cfg.solver exception=e
        end

        dt = time() - t0
        row = (
            cfg.name,
            String(cfg.solver),
            String(slab_type),
            String(slab_sizer),
            max_depth,
            warm_start_used,
            mip_seed_count,
            success,
            dt,
            n_beams,
            n_strips,
            norm_mass_beams,
            ec_total,
            max_deflection_in,
            err,
            beam_sizer_str,
            nlp_solver_out,
            deflection_limit_out,
            composite_action_out,
            staged_converged_out,
            staged_n_violations_out,
            n_L360_fail_out,
            n_L240_fail_out,
            global_δ_ok_out,
            max_staged_mm_out,
        )

        lock(table_lock) do
            push!(out_df, row)
            write_table_atomically(results_file, out_df)
        end
    end
end

function done_set_material_mc(results_file::String)::Set{Tuple{String, String, String}}
    out = Set{Tuple{String, String, String}}()
    if !isfile(results_file)
        return out
    end
    try
        df = CSV.read(results_file, DataFrame)
        required = [:geometry, :concrete_scenario, :slab_sizer]
        if !all(c -> c in names(df), required)
            return out
        end
        for r in eachrow(df)
            push!(out, (String(r.geometry), String(r.concrete_scenario), String(r.slab_sizer)))
        end
    catch e
        @warn "Failed to read material_scenario_mc done-set" file=results_file exception=e
    end
    return out
end

"""
    rand_lognormal(rng, nominal, σ_ln, N)

Draw N samples from a log-normal distribution with given nominal mean and
log-scale standard deviation. Avoids Distributions.jl dependency.
"""
function rand_lognormal(rng::AbstractRNG, nominal::Real, σ_ln::Real, N::Int)
    μ_ln = log(nominal) - σ_ln^2 / 2
    return exp.(μ_ln .+ σ_ln .* randn(rng, N))
end

"""
    mc_summary(samples) -> NamedTuple

Compute summary statistics for a 1-D vector of MC samples:
mean, std, p5/p25/p50/p75/p95 percentiles.
"""
function mc_summary(samples::AbstractVector{<:Real})
    return (
        mean = mean(samples),
        std  = std(samples),
        p5   = quantile(samples, 0.05),
        p25  = quantile(samples, 0.25),
        p50  = quantile(samples, 0.50),
        p75  = quantile(samples, 0.75),
        p95  = quantile(samples, 0.95),
    )
end

function run_material_scenario_mc(results_path::String; json_path::String=DEFAULT_TOPOLOGY_JSON_PATH)
    println("\n=== Study: material_scenario_mc ===")
    ensure_dir(results_path)

    results_file = joinpath(results_path, "material_scenario_mc.csv")
    done = done_set_material_mc(results_file)

    slab_sizers = [:uniform, :cellular]
    scenarios = collect(values(SlabDesignFactors.CONCRETE_SCENARIOS))
    sort!(scenarios, by=s -> s.name)

    configs = NamedTuple[]
    sub_paths = sort(filter(x -> endswith(x, ".json"), readdir(json_path)))
    for sub_path in sub_paths
        name = replace(sub_path, ".json" => "")
        path = joinpath(json_path, sub_path)
        for scenario in scenarios
            for sizer in slab_sizers
                key = (name, scenario.name, String(sizer))
                if key in done
                    continue
                end
                push!(configs, (name=name, path=path, scenario=scenario, slab_sizer=sizer))
            end
        end
    end

    println("Pending material_scenario_mc configs: $(length(configs))")
    isempty(configs) && return

    # Pre-compute MC log-scale sigmas
    σ_steel    = sqrt(log(1 + MC_CV_STEEL^2))
    σ_rebar    = sqrt(log(1 + MC_CV_REBAR^2))
    σ_concrete = sqrt(log(1 + MC_CV_CONCRETE^2))

    # Seed all columns so we can build the output DataFrame.
    out_df = if isfile(results_file)
        CSV.read(results_file, DataFrame)
    else
        DataFrame(
            geometry=String[], concrete_scenario=String[], slab_sizer=String[],
            success=Bool[], elapsed_s=Float64[],
            # structural
            ρ_concrete_kgm3=Float64[], E_c_ksi=Float64[], ecc_concrete_nom=Float64[],
            n_beams=Int[], n_strips=Int[],
            norm_mass_beams=Float64[], norm_mass_slab=Float64[], norm_mass_rebar=Float64[],
            ec_beams=Float64[], ec_slab=Float64[], ec_rebar=Float64[], ec_total=Float64[],
            max_deflection_in=Float64[],
            # MC — steel
            mc_steel_mean=Float64[], mc_steel_std=Float64[],
            mc_steel_p5=Float64[], mc_steel_p95=Float64[],
            # MC — concrete
            mc_conc_mean=Float64[], mc_conc_std=Float64[],
            mc_conc_p5=Float64[], mc_conc_p95=Float64[],
            # MC — rebar
            mc_rebar_mean=Float64[], mc_rebar_std=Float64[],
            mc_rebar_p5=Float64[], mc_rebar_p95=Float64[],
            # MC — joint (all three varied simultaneously)
            mc_joint_mean=Float64[], mc_joint_std=Float64[],
            mc_joint_p5=Float64[], mc_joint_p25=Float64[],
            mc_joint_p50=Float64[], mc_joint_p75=Float64[], mc_joint_p95=Float64[],
            # variance decomposition
            var_frac_steel=Float64[], var_frac_conc=Float64[], var_frac_rebar=Float64[],
            error_message=String[],
        )
    end
    ensure_results_columns!(out_df, STAGING_SUMMARY_COLS)

    table_lock = ReentrantLock()
    Threads.@threads for i in eachindex(configs)
        cfg = configs[i]
        t0 = time()
        scen = cfg.scenario
        println("[material_scenario_mc] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name) $(scen.name) $(cfg.slab_sizer)")

        success = false
        n_beams = 0; n_strips = 0
        norm_mass_beams = NaN; norm_mass_slab = NaN; norm_mass_rebar = NaN
        ec_beams = NaN; ec_slab = NaN; ec_rebar = NaN; ec_total = NaN
        max_δ = NaN
        mc_steel = (mean=NaN, std=NaN, p5=NaN, p25=NaN, p50=NaN, p75=NaN, p95=NaN)
        mc_conc  = (mean=NaN, std=NaN, p5=NaN, p25=NaN, p50=NaN, p75=NaN, p95=NaN)
        mc_reb   = (mean=NaN, std=NaN, p5=NaN, p25=NaN, p50=NaN, p75=NaN, p95=NaN)
        mc_joint = (mean=NaN, std=NaN, p5=NaN, p25=NaN, p50=NaN, p75=NaN, p95=NaN)
        vf_steel = NaN; vf_conc = NaN; vf_rebar = NaN
        err = ""

        beam_sizer_str = missing
        nlp_solver_out = missing
        deflection_limit_out = missing
        composite_action_out = missing
        staged_converged_out = missing
        staged_n_violations_out = missing
        n_L360_fail_out = missing
        n_L240_fail_out = missing
        global_δ_ok_out = missing
        max_staged_mm_out = missing

        try
            geometry_dict = read_geometry_from_json(cfg.path)
            geometry = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=false)

            slab_params = SlabDesignFactors.SlabAnalysisParams(
                geometry,
                slab_name=cfg.name,
                slab_type=DEFAULT_SLAB_TYPE,
                vector_1d=DEFAULT_VECTOR_1D,
                slab_sizer=cfg.slab_sizer,
                spacing=DEFAULT_SPACING,
                plot_analysis=false,
                fix_param=true,
                slab_units=:m,
            )

            sizing_params = SlabDesignFactors.SlabSizingParams(
                live_load=DEFAULT_LIVE_LOAD,
                superimposed_dead_load=DEFAULT_SDL,
                slab_dead_load=0.0,
                live_factor=DEFAULT_LIVE_FACTOR,
                dead_factor=DEFAULT_DEAD_FACTOR,
                beam_sizer=:discrete,
                max_depth=DEFAULT_MAX_DEPTH,
                beam_units=:in,
                serviceability_lim=DEFAULT_SERV_LIM,
                minimum_continuous=true,
                collinear=true,
                composite_action=true,
                deflection_reduction_factor=1.0,
                E_c=scen.E_c_ksi,
                concrete_material=scen,
            )

            slab_params = SlabDesignFactors.analyze_slab(slab_params)
            slab_params, sizing_params = SlabDesignFactors.optimal_beamsizer(slab_params, sizing_params)

            if !isempty(sizing_params.minimizers)
                results = SlabDesignFactors.postprocess_slab(slab_params, sizing_params, check_collinear=false)
                beam_sizer_str = String(results.beam_sizer)
                nlp_solver_out = results.nlp_solver
                deflection_limit_out = results.deflection_limit
                composite_action_out = results.composite_action
                staged_converged_out = results.staged_converged
                staged_n_violations_out = results.staged_n_violations
                n_L360_fail_out = results.n_L360_fail
                n_L240_fail_out = results.n_L240_fail
                global_δ_ok_out = results.global_δ_ok
                max_staged_mm_out = round(results.max_δ_total * 25.4, digits=3)
                n_beams = length(results.Δ_local)
                n_strips = length(slab_params.load_areas)
                norm_mass_beams = results.norm_mass_beams
                norm_mass_slab  = results.norm_mass_slab
                norm_mass_rebar = results.norm_mass_rebar
                ec_beams = results.embodied_carbon_beams
                ec_slab  = results.embodied_carbon_slab
                ec_rebar = results.embodied_carbon_rebar
                ec_total = ec_beams + ec_slab + ec_rebar
                max_δ = maximum(maximum(abs.(d)) for d in results.Δ_local if !isempty(d); init=0.0)

                # ── Monte Carlo ──────────────────────────────────────────
                # Per-config RNG (deterministic seed for reproducibility)
                rng = MersenneTwister(hash((cfg.name, scen.name, cfg.slab_sizer)))

                ecc_steel_samples    = rand_lognormal(rng, MC_ECC_STEEL_NOM,  σ_steel,    MC_N_SAMPLES)
                ecc_rebar_samples    = rand_lognormal(rng, MC_ECC_REBAR_NOM,  σ_rebar,    MC_N_SAMPLES)
                ecc_concrete_samples = rand_lognormal(rng, scen.ecc_concrete, σ_concrete, MC_N_SAMPLES)

                # Independent sweeps (only one coefficient varies per sweep)
                ec_steel_only    = norm_mass_beams .* ecc_steel_samples .+ ec_slab .+ ec_rebar
                ec_conc_only     = ec_beams .+ norm_mass_slab .* ecc_concrete_samples .+ ec_rebar
                ec_rebar_only    = ec_beams .+ ec_slab .+ norm_mass_rebar .* ecc_rebar_samples

                # Joint sweep (all three vary)
                ec_joint = norm_mass_beams .* ecc_steel_samples .+
                           norm_mass_slab  .* ecc_concrete_samples .+
                           norm_mass_rebar .* ecc_rebar_samples

                mc_steel = mc_summary(ec_steel_only)
                mc_conc  = mc_summary(ec_conc_only)
                mc_reb   = mc_summary(ec_rebar_only)
                mc_joint = mc_summary(ec_joint)

                # Variance decomposition (marginal variance contributions)
                v_steel = var(norm_mass_beams .* ecc_steel_samples)
                v_conc  = var(norm_mass_slab  .* ecc_concrete_samples)
                v_reb   = var(norm_mass_rebar .* ecc_rebar_samples)
                v_total = v_steel + v_conc + v_reb
                if v_total > 0
                    vf_steel = v_steel / v_total
                    vf_conc  = v_conc  / v_total
                    vf_rebar = v_reb   / v_total
                else
                    vf_steel = vf_conc = vf_rebar = 0.0
                end

                success = true
            else
                err = "No minimizers returned from optimizer"
            end
        catch e
            err = sprint(showerror, e)
            @warn "Failed material_scenario_mc config" name=cfg.name scenario=scen.name slab_sizer=cfg.slab_sizer exception=e
        end

        dt = time() - t0
        row = (
            cfg.name, scen.name, String(cfg.slab_sizer),
            success, dt,
            scen.ρ_concrete, scen.E_c_ksi, scen.ecc_concrete,
            n_beams, n_strips,
            norm_mass_beams, norm_mass_slab, norm_mass_rebar,
            ec_beams, ec_slab, ec_rebar, ec_total,
            max_δ,
            mc_steel.mean, mc_steel.std, mc_steel.p5, mc_steel.p95,
            mc_conc.mean,  mc_conc.std,  mc_conc.p5,  mc_conc.p95,
            mc_reb.mean,   mc_reb.std,   mc_reb.p5,   mc_reb.p95,
            mc_joint.mean, mc_joint.std, mc_joint.p5, mc_joint.p25,
            mc_joint.p50, mc_joint.p75, mc_joint.p95,
            vf_steel, vf_conc, vf_rebar,
            err,
            beam_sizer_str,
            nlp_solver_out,
            deflection_limit_out,
            composite_action_out,
            staged_converged_out,
            staged_n_violations_out,
            n_L360_fail_out,
            n_L240_fail_out,
            global_δ_ok_out,
            max_staged_mm_out,
        )

        lock(table_lock) do
            push!(out_df, row)
            write_table_atomically(results_file, out_df)
        end
    end
end

function done_set_validation_mip(results_file::String)::Set{String}
    out = Set{String}()
    if !isfile(results_file)
        return out
    end
    try
        df = CSV.read(results_file, DataFrame)
        if !("name" in names(df))
            return out
        end
        for r in eachrow(df)
            push!(out, String(r.name))
        end
    catch e
        @warn "Failed to read validation_mip done-set" file=results_file exception=e
    end
    return out
end

function run_validation_mip(results_path::String; json_path::String=DEFAULT_VALIDATION_JSON_PATH)
    println("\n=== Study: validation_mip ===")
    ensure_dir(results_path)

    results_name = "validation_mip"
    results_file = joinpath(results_path, results_name * ".csv")
    done = done_set_validation_mip(results_file)

    live_load_dict = Dict(
        "office" => 60,
        "school" => 50,
        "warehouse" => 250,
    )
    building_names = ["office", "school", "warehouse"]

    configs = NamedTuple[]
    for bname in building_names
        if bname in done
            continue
        end
        path = joinpath(json_path, bname * ".json")
        push!(configs, (name=bname, path=path))
    end

    println("Pending validation_mip configs: $(length(configs))")
    isempty(configs) && return

    append_lock = ReentrantLock()
    Threads.@threads for i in eachindex(configs)
        cfg = configs[i]
        println("[validation_mip] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name)")
        try
            geometry_dict = read_geometry_from_json(cfg.path)
            geometry, type_information = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=true)

            slab_params = SlabDesignFactors.SlabAnalysisParams(
                geometry,
                slab_name=cfg.name,
                slab_type=:uniaxial,
                slab_thickness=4.75 * SlabDesignFactors.convert_to_m[:in],
                vector_1d=[0.0, 1.0],
                slab_sizer=:uniform,
                spacing=DEFAULT_SPACING,
                plot_analysis=false,
                fix_param=true,
                slab_units=:m,
                i_holes=type_information["i_holes"],
                i_perimeter=type_information["i_perimeter"],
            )

            beam_sizing_params = SlabDesignFactors.SlabSizingParams(
                live_load=SlabDesignFactors.psf_to_ksi(live_load_dict[cfg.name]),
                superimposed_dead_load=SlabDesignFactors.psf_to_ksi(20),
                slab_dead_load=SlabDesignFactors.psf_to_ksi(45),
                façade_load=SlabDesignFactors.plf_to_kpi(500),
                live_factor=DEFAULT_LIVE_FACTOR,
                dead_factor=DEFAULT_DEAD_FACTOR,
                beam_sizer=:discrete,
                max_depth=DEFAULT_MAX_DEPTH,
                beam_units=:in,
                serviceability_lim=DEFAULT_SERV_LIM,
                collinear=false,
                drawn=true,
                element_ids=type_information["element_ids"],
                minimum_continuous=true,
                n_max_sections=0,
            )

            slab_params = SlabDesignFactors.analyze_slab(slab_params)
            slab_params, beam_sizing_params = SlabDesignFactors.optimal_beamsizer(slab_params, beam_sizing_params)
            if isempty(beam_sizing_params.minimizers)
                @warn "No minimizers for validation_mip" name=cfg.name
            else
                slab_results = SlabDesignFactors.postprocess_slab(slab_params, beam_sizing_params)
                lock(append_lock) do
                    SlabDesignFactors.append_results_to_csv(results_path * "/", results_name, [slab_results])
                end
            end
        catch e
            @warn "Failed validation_mip config" name=cfg.name exception=e
        end
    end
end

function analyze_experiments(
    results_path::String,
    completion_file::String;
    studies::Set{String}=DEFAULT_STUDIES,
)
    println("Dependencies loaded successfully.")
    println("Using $(Threads.nthreads()) Julia threads.")
    println("Selected studies: $(join(sort(collect(studies)), ", "))")
    BLAS.set_num_threads(1)

    ensure_dir(results_path)

    if "max_depths" in studies
        run_max_depths(results_path)
    end
    if "strip_resolution" in studies
        run_strip_resolution(results_path)
    end
    if "constrained_inventory" in studies
        run_constrained_inventory(results_path)
    end
    if "nlp_solver_comparison" in studies
        run_nlp_solver_comparison(results_path)
    end
    if "material_scenario_mc" in studies
        run_material_scenario_mc(results_path)
    end
    if "validation_mip" in studies
        run_validation_mip(results_path)
    end

    open(completion_file, "w") do f
        write(f, "Experiments complete")
    end
end

function main()
    args = ARGS
    if !(length(args) in (2, 3))
        println("Usage: julia executable_experiments.jl <results_path> <completion_file> [studies_csv]")
        println("studies_csv options: max_depths, strip_resolution, constrained_inventory, nlp_solver_comparison, material_scenario_mc, validation_mip")
        return
    end

    results_path = args[1]
    completion_file = args[2]
    studies = length(args) == 3 ? parse_studies_arg(args[3]) : DEFAULT_STUDIES

    analyze_experiments(results_path, completion_file; studies=studies)
end

main()
