using CSV
using DataFrames

"""
    discover_shard_csvs_nested(shards_dir) -> Dict{String, Vector{String}}

New layout: `shards/shard_NNN/<results_name>.csv`
"""
function discover_shard_csvs_nested(shards_dir::String)::Dict{String, Vector{String}}
    grouped = Dict{String, Vector{String}}()
    for name in sort(readdir(shards_dir))
        startswith(name, "shard_") || continue
        sub = joinpath(shards_dir, name)
        isdir(sub) || continue
        for file in readdir(sub)
            endswith(file, ".csv") || continue
            results_name = replace(file, ".csv" => "")
            push!(get!(grouped, results_name, String[]), joinpath(sub, file))
        end
    end
    return grouped
end

"""
    discover_shard_csvs_flat(shards_dir) -> Dict{String, Vector{String}}

Legacy layout: `shards/<results_name>_shard_NNN.csv`
"""
function discover_shard_csvs_flat(shards_dir::String)::Dict{String, Vector{String}}
    grouped = Dict{String, Vector{String}}()
    for file in sort(filter(f -> endswith(f, ".csv"), readdir(shards_dir)))
        m = match(r"^(.*)_shard_[0-9]+\.csv$", file)
        if m === nothing
            @warn "Skipping file that does not match shard pattern" file=file
            continue
        end
        results_name = m.captures[1]
        push!(get!(grouped, results_name, String[]), joinpath(shards_dir, file))
    end
    return grouped
end

function has_nested_shard_dirs(shards_dir::String)::Bool
    isdir(shards_dir) || return false
    return any(readdir(shards_dir)) do name
        startswith(name, "shard_") && isdir(joinpath(shards_dir, name))
    end
end

"""
    merge_shards(results_root, run_name, merged_dir)

Merge per-shard CSV outputs into one CSV per results name.

Supports:
- **Nested** (current): `shards/shard_001/grid.csv`, …
- **Flat** (legacy): `shards/grid_shard_001.csv`, …
"""
function merge_shards(results_root::String, run_name::String, merged_dir::String)
    shards_dir = joinpath(results_root, run_name, "shards")
    if !isdir(shards_dir)
        error("Shards directory not found: $shards_dir")
    end
    mkpath(merged_dir)

    grouped = if has_nested_shard_dirs(shards_dir)
        discover_shard_csvs_nested(shards_dir)
    else
        discover_shard_csvs_flat(shards_dir)
    end

    if isempty(grouped)
        println("No shard CSV files found in $shards_dir")
        return
    end

    dedupe_cols = [
        :name, :slab_type, :slab_sizer, :beam_sizer, :collinear,
        :vector_1d_x, :vector_1d_y, :max_depth, :unique_sections,
    ]

    for (results_name, files) in sort(collect(grouped))
        isempty(files) && continue
        dfs = DataFrame[]
        for file in sort(files)
            push!(dfs, CSV.read(file, DataFrame))
        end
        merged = vcat(dfs..., cols=:union)

        active_dedupe_cols = filter(c -> c in names(merged), dedupe_cols)
        if !isempty(active_dedupe_cols)
            merged = unique(merged, active_dedupe_cols)
        else
            merged = unique(merged)
        end

        out_path = joinpath(merged_dir, results_name * ".csv")
        CSV.write(out_path, merged)
        println("Wrote merged file: $out_path (rows=$(nrow(merged)))")
    end
end

function main()
    args = ARGS
    if length(args) != 3
        println("Usage: julia executable_merge_shards.jl <results_root> <run_name> <merged_dir>")
        return
    end
    merge_shards(args[1], args[2], args[3])
end

main()
