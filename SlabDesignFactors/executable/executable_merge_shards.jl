using CSV
using DataFrames

"""
    merge_shards(results_root, run_name, merged_dir)

Merge per-shard CSV outputs into one CSV per results name.
"""
function merge_shards(results_root::String, run_name::String, merged_dir::String)
    shards_dir = joinpath(results_root, run_name, "shards")
    if !isdir(shards_dir)
        error("Shards directory not found: $shards_dir")
    end
    mkpath(merged_dir)

    shard_files = sort(filter(f -> endswith(f, ".csv"), readdir(shards_dir)))
    if isempty(shard_files)
        println("No shard CSV files found in $shards_dir")
        return
    end

    grouped = Dict{String, Vector{String}}()
    for file in shard_files
        m = match(r"^(.*)_shard_[0-9]+\.csv$", file)
        if m === nothing
            @warn "Skipping file that does not match shard pattern" file=file
            continue
        end
        results_name = m.captures[1]
        get!(grouped, results_name, String[])
        push!(grouped[results_name], joinpath(shards_dir, file))
    end

    dedupe_cols = [
        :name, :slab_type, :slab_sizer, :beam_sizer, :collinear,
        :vector_1d_x, :vector_1d_y, :max_depth, :unique_sections,
    ]

    for (results_name, files) in sort(collect(grouped))
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
