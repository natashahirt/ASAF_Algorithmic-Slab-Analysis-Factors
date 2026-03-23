include("../SlabDesignFactors.jl")

using Base.Threads
using LinearAlgebra
using JSON
using CSV
using DataFrames
using .SlabDesignFactors

"""
    parse_params_file(params_file::String) -> Vector{Tuple{String, String}}

Read executable parameter lines of the form:
`<json_path> <results_name>`.
Blank lines and `#` comments are ignored.
"""
function parse_params_file(params_file::String)::Vector{Tuple{String, String}}
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
        push!(entries, (parts[1], parts[2]))
    end
    return entries
end

"""
    config_key(name, slab_sizer, max_depth, slab_type, vector_1d)

Canonical config identity used for resume logic.
"""
function config_key(
    name::String,
    slab_sizer::Symbol,
    max_depth::Int,
    slab_type::Symbol,
    vector_1d::Vector{Float64},
)::NTuple{6, Any}
    return (
        name,
        String(slab_sizer),
        Float64(max_depth),
        String(slab_type),
        Float64(vector_1d[1]),
        Float64(vector_1d[2]),
    )
end

"""
    done_set_for_file(results_file::String) -> Set

Load already-completed config keys from an existing shard CSV.
If the file is missing or unreadable, returns an empty set.
"""
function done_set_for_file(results_file::String)::Set{NTuple{6, Any}}
    out = Set{NTuple{6, Any}}()
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

    configs = NamedTuple[]
    for (json_path, results_name) in params_entries
        sub_paths = sort(filter(x -> endswith(x, ".json"), readdir(json_path)))
        for sub_path in sub_paths
            path = joinpath(json_path, sub_path)
            name = replace(sub_path, ".json" => "")
            for max_depth in max_depths
                for slab_sizer in slab_sizers
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
                        ))
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
    run_one_config(cfg) -> Union{Vector{SlabOptimResults}, Nothing}

Run a single slab configuration. Returns `nothing` on failure.
"""
function run_one_config(cfg)::Union{Vector{SlabOptimResults}, Nothing}
    try
        json_string = replace(read(cfg.path, String), "\\n" => "")
        geometry_dict = JSON.parse(JSON.parse(json_string, dicttype=Dict))
        geometry = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=false)

        slab_params = SlabDesignFactors.SlabAnalysisParams(
            geometry,
            slab_name=cfg.name,
            slab_type=cfg.slab_type,
            vector_1d=cfg.vector_1d,
            slab_sizer=cfg.slab_sizer,
            spacing=0.1,
            plot_analysis=false,
            fix_param=true,
            slab_units=:m,
        )

        beam_sizing_params = SlabDesignFactors.SlabSizingParams(
            live_load=SlabDesignFactors.psf_to_ksi(50),
            superimposed_dead_load=SlabDesignFactors.psf_to_ksi(15),
            live_factor=1.6,
            dead_factor=1.2,
            beam_sizer=:discrete,
            max_depth=cfg.max_depth,
            beam_units=:in,
            serviceability_lim=360,
            minimum_continuous=true,
            deflection_limit=false,
        )

        return collect(SlabDesignFactors.iterate_discrete_continuous(slab_params, beam_sizing_params))
    catch e
        @warn "Config failed; skipping" name=cfg.name slab_type=cfg.slab_type slab_sizer=cfg.slab_sizer max_depth=cfg.max_depth exception=e
        return nothing
    end
end

"""
    run_shard(results_root, run_name, completion_dir, params_file)

Run a deterministic shard of the full sweep, resumable from per-shard CSVs.
Shard metadata is read from Slurm environment variables.
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
    progress_dir = joinpath(results_root, run_name, "progress")
    completion_shards_dir = joinpath(completion_dir, run_name)
    mkpath(results_dir)
    mkpath(progress_dir)
    mkpath(completion_shards_dir)

    params_entries = parse_params_file(params_file)
    all_configs = build_all_configs(params_entries)
    shard_configs = [cfg for (i, cfg) in enumerate(all_configs) if mod(i - 1, shard_count) == shard_id - 1]

    # Resume state is isolated per {results_name, shard}.
    done_by_result = Dict{String, Set{NTuple{6, Any}}}()
    for (_, results_name) in params_entries
        shard_results_file = joinpath(results_dir, "$(results_name)_shard_$(shard_suffix).csv")
        done_by_result[results_name] = done_set_for_file(shard_results_file)
    end

    pending = NamedTuple[]
    for cfg in shard_configs
        key = config_key(cfg.name, cfg.slab_sizer, cfg.max_depth, cfg.slab_type, cfg.vector_1d)
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
        println("[Shard $(shard_id)] ($(Threads.threadid())/$(Threads.nthreads())) $(cfg.name) $(cfg.slab_type) $(cfg.vector_1d) $(cfg.slab_sizer) $(cfg.max_depth)in")

        iteration_result = run_one_config(cfg)
        if iteration_result === nothing
            Threads.atomic_add!(failed, 1)
            continue
        end

        lock(append_lock) do
            shard_filename = "$(cfg.results_name)_shard_$(shard_suffix)"
            SlabDesignFactors.append_results_to_csv(results_dir * "/", shard_filename, iteration_result)
            key = config_key(cfg.name, cfg.slab_sizer, cfg.max_depth, cfg.slab_type, cfg.vector_1d)
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
    params_file = length(args) == 4 ? args[4] : "/home/nhirt/2024_Slab-Design-Factors/SlabDesignFactors/executable/params.txt"

    run_shard(results_root, run_name, completion_dir, params_file)
end

main()
