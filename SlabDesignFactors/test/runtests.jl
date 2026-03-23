"""
Master test runner for SlabDesignFactors.

Usage (from project root):
    julia --project=. SlabDesignFactors/test/run.jl

Two-tier test suite:
  1. Unit tests   — fast, no geometry loading
  2. Integration  — full pipeline on real topology JSONs (slower)
"""

using Test

@testset "SlabDesignFactors" begin
    @testset "Unit tests" begin
        include("test_structs.jl")
        include("test_compression.jl")
        include("test_staged_deflection.jl")
    end
    @testset "Integration tests" begin
        include("test_integration.jl")
    end
end
