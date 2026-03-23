"""
Unit tests for SlabOptimResults @kwdef struct and related utilities.

Run from project root after loading the full module:
    include("SlabDesignFactors/SlabDesignFactors.jl")
    include("SlabDesignFactors/test/test_structs.jl")

These tests verify the @kwdef refactor works correctly: default construction,
keyword construction, field access, and the new δ_beam_dead field.
"""

using Test

@testset "SlabOptimResults @kwdef" begin

    @testset "default construction" begin
        r = SlabOptimResults()
        @test r.slab_name == ""
        @test r.slab_type == :uniaxial
        @test r.vector_1d == [1.0, 0.0]
        @test r.area == 0.0
        @test r.max_depth == 0.0
        @test r.collinear === nothing
        @test r.mass_beams == 0.0
        @test r.mass_columns == 0.0
        @test isempty(r.col_sections)
        @test isempty(r.δ_slab_dead)
        @test isempty(r.δ_beam_dead)
        @test isempty(r.δ_live_ok)
    end

    @testset "keyword construction" begin
        r = SlabOptimResults(
            slab_name   = "test_slab",
            area        = 100.0,
            max_depth   = 24.0,
            collinear   = true,
            beam_sizer  = :discrete,
            mass_beams  = 500.0,
            mass_columns= 200.0,
        )
        @test r.slab_name == "test_slab"
        @test r.area == 100.0
        @test r.max_depth == 24.0
        @test r.collinear == true
        @test r.beam_sizer == :discrete
        @test r.mass_beams == 500.0
        @test r.mass_columns == 200.0
        @test r.slab_type == :uniaxial
        @test r.embodied_carbon_beams == 0.0
        @test isempty(r.δ_beam_dead)
    end

    @testset "mutability" begin
        r = SlabOptimResults()
        r.slab_name = "mutated"
        r.δ_beam_dead = [0.001, 0.002, 0.003]
        @test r.slab_name == "mutated"
        @test length(r.δ_beam_dead) == 3
        @test r.δ_beam_dead[2] ≈ 0.002
    end

    @testset "staged deflection fields" begin
        n = 4
        r = SlabOptimResults(
            δ_slab_dead  = fill(0.005, n),
            δ_beam_dead  = fill(0.001, n),
            δ_sdl        = fill(0.003, n),
            δ_live       = fill(0.004, n),
            δ_total      = fill(0.013, n),
            Δ_limit_live = fill(0.02, n),
            Δ_limit_total= fill(0.03, n),
            δ_live_ok    = fill(true, n),
            δ_total_ok   = fill(true, n),
        )
        @test length(r.δ_slab_dead) == n
        @test length(r.δ_beam_dead) == n
        @test all(r.δ_live_ok)
        @test r.δ_total[1] ≈ 0.013
    end

    @testset "column sizing fields" begin
        r = SlabOptimResults(
            col_sections           = ["W10x49", "W12x53"],
            col_Pu                 = [100.0, 150.0],
            col_ϕPn                = [200.0, 250.0],
            col_util               = [0.5, 0.6],
            mass_columns           = 800.0,
            norm_mass_columns      = 5.0,
            embodied_carbon_columns= 6.1,
        )
        @test length(r.col_sections) == 2
        @test r.col_sections[1] == "W10x49"
        @test r.col_util[2] ≈ 0.6
        @test r.embodied_carbon_columns ≈ 6.1
    end
end

println("All SlabOptimResults tests passed.")
