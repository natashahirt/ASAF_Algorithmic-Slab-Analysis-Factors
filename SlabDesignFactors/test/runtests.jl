"""
Master test runner for SlabDesignFactors.

Usage (from project root):
    julia --project=. SlabDesignFactors/test/run.jl

Test tiers:
  1. Unit tests (`test_*.jl` in the unit block) — fast
  2. Integration — MIP/NLP sizing on topology JSONs (`test_integration.jl`)
  3. Slurm parity — `test_slurm_pipeline.jl` (`run_shard` + merge); uses `secrets/gurobi.lic` when present
"""

include(joinpath(@__DIR__, "_license_env.jl"))

using Test

@testset "SlabDesignFactors" begin
    @testset "Unit tests" begin
        include("test_structs.jl")
        include("test_compression.jl")
        include("test_slab_depth_minimum.jl")
        include("test_staged_deflection.jl")
        include("test_consolidate_loads.jl")
        include("test_results_csv.jl")
    end
    @testset "Integration tests" begin
        include("test_integration.jl")
        include("test_slurm_pipeline.jl")
    end
end
