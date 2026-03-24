"""
Cluster-parity smoke tests: same code paths as Slurm `executable_analyze_sharded.jl` + merge.

- Loads the sharded driver via `include` (CLI `main()` is skipped when not the entry script).
- Verifies double-encoded geometry JSON handling (`geometry_dict_from_json_path`).
- With a Gurobi license, runs `run_shard` on one topology JSON and `merge_shards`.

Run via `julia --project=. SlabDesignFactors/test/run.jl` (uses `_scripts.jl` first).
`run.jl` / `runtests.jl` set `GRB_LICENSE_FILE` to `secrets/gurobi.lic` when that file exists.
"""

using Test
using Gurobi

function _gurobi_license_available()
    try
        Gurobi.Env()
        return true
    catch
        return false
    end
end

const _SLURM_PIPELINE_GUROBI = _gurobi_license_available()
const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "Slurm pipeline smoke" begin
    @testset "Double-encoded geometry JSON → Dict & Asap.Model" begin
        grid_json = joinpath(_REPO_ROOT, "Geometries", "grid", "x0y0.json")
        @test isfile(grid_json)
        d = geometry_dict_from_json_path(grid_json)
        @test d isa Dict
        geom, meta = generate_from_json(d; plot=false, drawn=false)
        @test geom isa Asap.Model
        @test meta isa Dict && haskey(meta, "i_holes")
    end

    @testset "Include cluster drivers (Slurm code path; no CLI main())" begin
        # Same files as batch jobs; `main()` runs only when PROGRAM_FILE` is this script.
        include(joinpath(@__DIR__, "..", "executable", "executable_analyze_sharded.jl"))
        include(joinpath(@__DIR__, "..", "executable", "executable_merge_shards.jl"))
        @test isdefined(Main, :run_shard)
        @test isdefined(Main, :run_one_config)
        @test isdefined(Main, :merge_shards)
    end

    if !_SLURM_PIPELINE_GUROBI
        @info "Slurm pipeline: skipping run_shard + merge_shards (no Gurobi license). Add `secrets/gurobi.lic` or set GRB_LICENSE_FILE."
    else
        @testset "run_shard + merge_shards (one JSON, SLURM_ARRAY_TASK_COUNT=1)" begin
            mktempdir() do tmp
                geomdir = joinpath(tmp, "geometry_single")
                mkpath(geomdir)
                cp(
                    joinpath(_REPO_ROOT, "Geometries", "topology", "r1c1.json"),
                    joinpath(geomdir, "r1c1.json"),
                )
                params_file = joinpath(tmp, "params.txt")
                write(params_file, "$(geomdir) topology\n")

                results_root = joinpath(tmp, "results")
                completion_dir = joinpath(tmp, "completion")
                run_name = "smoke_run"
                merged_dir = joinpath(tmp, "merged")

                withenv(
                    "SLURM_ARRAY_TASK_ID" => "1",
                    "SLURM_ARRAY_TASK_COUNT" => "1",
                    "SLABDESIGN_SKIP_MAKIE" => "1",
                ) do
                    run_shard(results_root, run_name, completion_dir, params_file)
                end

                shard_mark = joinpath(completion_dir, run_name, "shard_001.complete")
                topology_csv = joinpath(results_root, run_name, "shards", "shard_001", "topology.csv")
                @test isfile(shard_mark)
                @test isfile(topology_csv)

                merge_shards(results_root, run_name, merged_dir)
                merged_csv = joinpath(merged_dir, "topology.csv")
                @test isfile(merged_csv)
            end
        end
    end
end
