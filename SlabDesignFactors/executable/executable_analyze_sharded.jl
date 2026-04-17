include("../SlabDesignFactors.jl")

using Base.Threads
using LinearAlgebra
using JSON
using Tables
using CSV
using DataFrames
using .SlabDesignFactors

"""
    parse_params_file(params_file::String) -> Vector{Tuple{String, String}}

Read lines `<geometry_dir> <results_name>`. Blank lines and `#` comments are ignored.

If `geometry_dir` is not an absolute path, it is resolved relative to the **repository root**:
`normpath(joinpath(dirname(abspath(params_file)), "..", ".."))` when `params_file` is under
`SlabDesignFactors/executable/`.
"""
function parse_params_file(params_file::String)::Vector{Tuple{String, String}}
    params_abs = abspath(params_file)
    repo_root = normpath(joinpath(dirname(params_abs), "..", ".."))
    entries = Tuple{String, String}[]
    for raw_line in readlines(params_file)
        line = strip(raw_line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        parts = split(line)
        if length(parts) < 2
            @warn "Skipping malformed params line" line=line
            continue
        end
        raw_path = String(parts[1])
        geom_dir = isabspath(raw_path) ? normpath(raw_path) : normpath(joinpath(repo_root, raw_path))
        push!(entries, (geom_dir, String(parts[2])))
    end
    return entries
end

"""
    config_key(name, slab_sizer, max_depth, slab_type, vector_1d, collinear, config_hash)

Canonical config identity used for resume logic.  The last element is the
`FULL_SWEEP_CONFIG_HASH` so that stale results produced by an older code
version or different parameter set are automatically re-run.
"""
function config_key(
    name::String,
    slab_sizer::Symbol,
    max_depth::Real,
    slab_type::Symbol,
    vector_1d::Vector{Float64},
    collinear::Bool,
    config_hash::String=SlabDesignFactors.FULL_SWEEP_CONFIG_HASH,
)::NTuple{8, Any}
    return (
        name,
        String(slab_sizer),
        Float64(max_depth),
        String(slab_type),
        Float64(vector_1d[1]),
        Float64(vector_1d[2]),
        collinear,
        config_hash,
    )
end

"""
    done_set_for_file(results_file::String) -> Set

Load already-completed config keys from an existing shard CSV.
If the file is missing or unreadable, returns an empty set.

When the CSV lacks a `config_hash` column (legacy data), the row is
treated as hash-less and will never match the current hash, forcing a
re-run — which is exactly the desired behaviour after a code change.
"""
function done_set_for_file(results_file::String)::Set{NTuple{8, Any}}
    out = Set{NTuple{8, Any}}()
    if !isfile(results_file)
        return out
    end
    try
        df = CSV.read(results_file, DataFrame)
        required_cols = (:name, :slab_sizer, :max_depth, :slab_type, :vector_1d_x, :vector_1d_y)
        if !all(col -> col in names(df), required_cols)
            @warn "Resume CSV is missing expected columns; continuing with empty done-set" file=results_file
            return out
        end
        if !("collinear" in names(df))
            @warn "Resume CSV lacks collinear column; clearing done-set so both collinearity sweeps re-run" file=results_file
            return out
        end
        has_hash = "config_hash" in names(df)
        for r in eachrow(df)
            row_hash = has_hash ? String(r.config_hash) : ""
            push!(out, (
                String(r.name),
                String(r.slab_sizer),
                Float64(r.max_depth),
                String(r.slab_type),
                Float64(r.vector_1d_x),
                Float64(r.vector_1d_y),
                Bool(r.collinear),
                row_hash,
            ))
        end
    catch e
        @warn "Failed to read done-set CSV; continuing from empty set" file=results_file exception=e
    end
    return out
end

"""
    build_all_configs(params_entries) -> Vector{NamedTuple}

Generate the full sweep config list across all geometry groups.
"""
function build_all_configs(params_entries::Vector{Tuple{String, String}})
    slab_types = [:isotropic, :orth_biaxial, :orth_biaxial, :uniaxial, :uniaxial, :uniaxial, :uniaxial]
    vector_1ds = [
        [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, -1.0],
    ]
    max_depths = [25, 40]
    slab_sizers = [:cellular, :uniform]
    collinears = [true, false]

    configs = NamedTuple[]
    for (json_path, results_name) in params_entries
        sub_paths = sort(filter(x -> endswith(x, ".json"), readdir(json_path)))
        for sub_path in sub_paths
            path = joinpath(json_path, sub_path)
            name = replace(sub_path, ".json" => "")
            for max_depth in max_depths
                for slab_sizer in slab_sizers
                    for collinear in collinears
                        for (i, slab_type) in enumerate(slab_types)
                            vector_1d = vector_1ds[i]
                            push!(configs, (
                                json_path=json_path,
                                results_name=results_name,
                                path=path,
                                name=name,
                                slab_type=slab_type,
                                vector_1d=vector_1d,
                                slab_sizer=slab_sizer,
                                max_depth=max_depth,
                                collinear=collinear,
                            ))
                        end
                    end
                end
            end
        end
    end
    return configs
end

"""
    write_progress(progress_file, payload)

Write a JSON progress snapshot atomically.
"""
function write_progress(progress_file::String, payload::Dict{String, Any})
    tmp_file = progress_file * ".tmp"
    open(tmp_file, "w") do io
        JSON.print(io, payload)
    end
    Base.Filesystem.rename(tmp_file, progress_file)
end

"""
    run_one_config(cfg) -> Union{Vector{SlabDesignFactors.SlabOptimResults}, Nothing}

Run a single slab configuration. Returns `nothing` on failure.
"""
function run_one_config(cfg)::Union{Vector{SlabDesignFactors.SlabOptimResults}, Nothing}
    try
        geometry_dict = SlabDesignFactors.geometry_dict_from_json_path(cfg.path)
        geometry, _ = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=false)

        slab_params = SlabDesignFactors.SlabAnalysisParams(
            geometry,
            slab_name=cfg.name,
            slab_type=cfg.slab_type,
            vector_1d=cfg.vector_1d,
            slab_sizer=cfg.slab_sizer,
            spacing=SlabDesignFactors.FULL_SWEEP_STRIP_SPACING,
            plot_analysis=false,
            fix_param=true,
            slab_units=:m,
        )

        beam_sizing_params = SlabDesignFactors.SlabSizingParams(
            live_load=SlabDesignFactors.FULL_SWEEP_LIVE_LOAD,
            superimposed_dead_load=SlabDesignFactors.FULL_SWEEP_SUPERIMPOSED_DEAD_LOAD,
            live_factor=SlabDesignFactors.FULL_SWEEP_LIVE_FACTOR,
            dead_factor=SlabDesignFactors.FULL_SWEEP_DEAD_FACTOR,
            beam_sizer=:discrete,
            max_depth=cfg.max_depth,
            beam_units=:in,
            serviceability_lim=SlabDesignFactors.FULL_SWEEP_SERVICEABILITY_LIM,
            minimum_continuous=true,
            collinear=cfg.collinear,
            composite_action=SlabDesignFactors.FULL_SWEEP_COMPOSITE_ACTION,
            staged_deflection_limit=SlabDesignFactors.FULL_SWEEP_STAGED_DEFLECTION_LIMIT,
            deflection_reduction_factor=SlabDesignFactors.FULL_SWEEP_DEFLECTION_REDUCTION_FACTOR,
        )

        results = collect(SlabDesignFactors.iterate_discrete_continuous(slab_params, beam_sizing_params))
        for r in results
            r.config_hash = SlabDesignFactors.FULL_SWEEP_CONFIG_HASH
        end
        return results
    catch e
        @warn "Config failed; skipping" name=cfg.name slab_type=cfg.slab_type slab_sizer=cfg.slab_sizer max_depth=cfg.max_depth collinear=cfg.collinear exception=e
        return nothing
    end
end

"""
    run_shard(results_root, run_name, completion_dir, params_file)

Run a deterministic shard of the full sweep, resumable from per-shard CSVs.
Shard metadata is read from Slurm environment variables.

Outputs for this shard go under `results_root/run_name/shards/shard_<NNN>/` (e.g. `grid.csv`, `nova.csv`).
Progress JSON stays in `results_root/run_name/progress/`.
"""
function run_shard(
    results_root::String,
    run_name::String,
    completion_dir::String,
    params_file::String,
)
    shard_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
    shard_count = parse(Int, get(ENV, "SLURM_ARRAY_TASK_COUNT", "1"))
    shard_suffix = lpad(string(shard_id), 3, '0')

    BLAS.set_num_threads(1)
    println("Shard $(shard_id)/$(shard_count) starting with $(Threads.nthreads()) Julia threads and 1 BLAS thread.")

    results_dir = joinpath(results_root, run_name, "shards")
    shard_dir = joinpath(results_dir, "shard_$(shard_suffix)")
    progress_dir = joinpath(results_root, run_name, "progress")
    completion_shards_dir = joinpath(completion_dir, run_name)
    mkpath(shard_dir)
    mkpath(progress_dir)
    mkpath(completion_shards_dir)

    params_entries = parse_params_file(params_file)
    all_configs = build_all_configs(params_entries)
    shard_configs = [cfg for (i, cfg) in enumerate(all_configs) if mod(i - 1, shard_count) == shard_id - 1]

    # Resume state is isolated per {results_name, shard}.
    done_by_result = Dict{String, Set{NTuple{8, Any}}}()
    for (_, results_name) in params_entries
        shard_results_file = joinpath(shard_dir, "$(results_name).csv")
        done_by_result[results_name] = done_set_for_file(shard_results_file)
    end

    pending = NamedTuple[]
    for cfg in shard_configs
        key = config_key(cfg.name, cfg.slab_sizer, cfg.max_depth, cfg.slab_type, cfg.vector_1d, cfg.collinear)
        if !(key in done_by_result[cfg.results_name])
            push!(pending, cfg)
        end
    end

    println("Total configs: $(length(all_configs)) | shard configs: $(length(shard_configs)) | pending: $(length(pending))")

    append_lock = ReentrantLock()
    progress_lock = ReentrantLock()
    completed = Threads.Atomic{Int}(0)
    failed = Threads.Atomic{Int}(0)

    progress_file = joinpath(progress_dir, "progress_shard_$(shard_suffix).json")
    write_progress(progress_file, Dict(
        "shard_id" => shard_id,
        "shard_count" => shard_count,
        "completed" => 0,
        "failed" => 0,
        "pending" => length(pending),
        "total_shard_configs" => length(shard_configs),
        "status" => "running",
    ))

    Threads.@threads for idx in eachindex(pending)
        cfg = pending[idx]
        println("[Shard $(shard_id)] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name) $(cfg.slab_type) $(cfg.vector_1d) $(cfg.slab_sizer) $(cfg.max_depth)in collinear=$(cfg.collinear)")

        iteration_result = run_one_config(cfg)
        if iteration_result === nothing
            Threads.atomic_add!(failed, 1)
            continue
        end

        lock(append_lock) do
            SlabDesignFactors.append_results_to_csv(shard_dir * "/", cfg.results_name, iteration_result)
            key = config_key(cfg.name, cfg.slab_sizer, cfg.max_depth, cfg.slab_type, cfg.vector_1d, cfg.collinear)
            push!(done_by_result[cfg.results_name], key)
        end

        new_completed = Threads.atomic_add!(completed, 1) + 1
        if new_completed % 25 == 0 || new_completed == length(pending)
            lock(progress_lock) do
                write_progress(progress_file, Dict(
                    "shard_id" => shard_id,
                    "shard_count" => shard_count,
                    "completed" => new_completed,
                    "failed" => failed[],
                    "pending" => max(length(pending) - new_completed, 0),
                    "total_shard_configs" => length(shard_configs),
                    "status" => "running",
                ))
            end
        end
    end

    write_progress(progress_file, Dict(
        "shard_id" => shard_id,
        "shard_count" => shard_count,
        "completed" => completed[],
        "failed" => failed[],
        "pending" => max(length(pending) - completed[], 0),
        "total_shard_configs" => length(shard_configs),
        "status" => "finished",
    ))

    shard_completion_file = joinpath(completion_shards_dir, "shard_$(shard_suffix).complete")
    open(shard_completion_file, "w") do io
        write(io, "completed")
    end

    markers = filter(x -> startswith(x, "shard_") && endswith(x, ".complete"), readdir(completion_shards_dir))
    if length(markers) == shard_count
        run_completion_file = joinpath(completion_shards_dir, "run_complete.txt")
        open(run_completion_file, "w") do io
            write(io, "all shards completed")
        end
    end

    println("Shard $(shard_id) done. Completed=$(completed[]) Failed=$(failed[])")
end

function main()
    args = ARGS
    if !(length(args) in (3, 4))
        println("Usage: julia executable_analyze_sharded.jl <results_root> <run_name> <completion_dir> [params_file]")
        return
    end

    results_root = args[1]
    run_name = args[2]
    completion_dir = args[3]
    params_file = length(args) == 4 ? args[4] : joinpath(@__DIR__, "params.txt")

    run_shard(results_root, run_name, completion_dir, params_file)
end

# Run CLI only when this file is the Julia entry script (not when `include`d from tests).
if !isempty(PROGRAM_FILE) && abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
