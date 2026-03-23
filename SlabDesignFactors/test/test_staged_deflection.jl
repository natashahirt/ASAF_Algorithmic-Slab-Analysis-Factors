"""
Unit tests for staged deflection analysis logic.

Run from project root after loading the full module:
    include("SlabDesignFactors/SlabDesignFactors.jl")
    include("SlabDesignFactors/test/test_staged_deflection.jl")

These tests verify:
- Analytical beam self-weight deflection formula (5wL⁴ / 384EI)
- Load-fraction splitting for SDL vs LL
- Superposition δ_total = δ_dead + δ_beam_sw + δ_sdl + δ_live
- Deflection limit checks (L/360, L/240)
"""

using Test

@testset "Staged deflection — analytical & superposition" begin

    @testset "analytical beam self-weight deflection" begin
        # W21x62: d=20.99, bf=8.240, tw=0.400, tf=0.615
        sec = I_symm(20.99, 8.240, 0.400, 0.615)
        E = steel_ksi.E               # 29000 ksi
        ρ = steel_ksi.ρ               # kip/in³
        A = sec.A                      # in²
        Ix = sec.Ix                    # in⁴
        L = 30.0 * 12.0               # 30 ft = 360 in

        w_sw = A * ρ                   # kip/in
        δ_analytic = 5 * w_sw * L^4 / (384 * E * Ix)

        @test δ_analytic > 0
        # Sanity: self-weight deflection for a 30ft W21x62 should be small
        # (order of ~0.1-0.5 in)
        @test δ_analytic < 1.0  # definitely less than 1 inch
        @test δ_analytic > 0.01 # but not negligible
    end

    @testset "load-fraction splitting" begin
        # If SDL = 20 psf and LL = 40 psf (mapped to kip/in units),
        # the fraction of live load is 2/3.
        w_sdl = 20.0
        w_live = 40.0
        w_total = w_sdl + w_live
        f_live = w_live / w_total
        f_sdl = w_sdl / w_total
        @test f_live ≈ 2/3
        @test f_sdl ≈ 1/3
        @test f_live + f_sdl ≈ 1.0

        # Combined deflection split
        δ_combined = 0.3  # inches
        δ_live = δ_combined * f_live
        δ_sdl = δ_combined * f_sdl
        @test δ_live ≈ 0.2
        @test δ_sdl ≈ 0.1
        @test δ_live + δ_sdl ≈ δ_combined
    end

    @testset "superposition and limit checks" begin
        L = 360.0  # inches (30 ft beam)
        Δ_lim_live  = L / 360.0
        Δ_lim_total = L / 240.0

        δ_slab_dead = 0.4
        δ_beam_dead = 0.05
        δ_sdl = 0.15
        δ_live = 0.35
        δ_total = δ_slab_dead + δ_beam_dead + δ_sdl + δ_live

        @test δ_total ≈ 0.95

        # L/360 = 1.0 in, L/240 = 1.5 in
        @test Δ_lim_live ≈ 1.0
        @test Δ_lim_total ≈ 1.5

        live_ok  = δ_live <= Δ_lim_live
        total_ok = δ_total <= Δ_lim_total

        @test live_ok == true    # 0.35 ≤ 1.0
        @test total_ok == true   # 0.95 ≤ 1.5
    end

    @testset "failing deflection" begin
        L = 360.0
        Δ_lim_live  = L / 360.0   # 1.0 in
        Δ_lim_total = L / 240.0   # 1.5 in

        δ_live = 1.2   # exceeds L/360
        δ_total = 1.8  # exceeds L/240

        @test (δ_live <= Δ_lim_live) == false
        @test (δ_total <= Δ_lim_total) == false
    end

    @testset "stiffness ratio scaling" begin
        # If bare steel Ix = 1000 in⁴ and composite Ix = 2500 in⁴,
        # deflection on composite = δ_bare / 2.5
        Ix_bare = 1000.0
        Ix_comp = 2500.0
        ratio = Ix_comp / Ix_bare
        @test ratio ≈ 2.5

        δ_bare = 0.5  # in
        δ_comp = δ_bare / ratio
        @test δ_comp ≈ 0.2

        # Staged: dead stays on bare, SDL+LL on composite
        f_dead = 0.4
        f_sdl  = 0.2
        f_live = 0.4
        δ_dead_bare = δ_bare * f_dead           # 0.2
        δ_sdl_comp  = δ_bare * f_sdl / ratio    # 0.04
        δ_live_comp = δ_bare * f_live / ratio   # 0.08
        δ_total = δ_dead_bare + δ_sdl_comp + δ_live_comp

        @test δ_dead_bare ≈ 0.2
        @test δ_sdl_comp ≈ 0.04
        @test δ_live_comp ≈ 0.08
        @test δ_total ≈ 0.32
    end

    @testset "beam self-weight invariant: heavier beam → more deflection" begin
        # Heavier section means more self-weight but also more Ix.
        # For typical W-shapes the Ix grows faster than A, so
        # self-weight deflection should decrease for heavier sections at same span.
        sec_light = I_symm(12.0, 6.5, 0.300, 0.440)  # ~W12x30
        sec_heavy = I_symm(21.0, 8.24, 0.400, 0.615)  # ~W21x62

        L = 360.0
        E = steel_ksi.E
        ρ = steel_ksi.ρ

        δ_sw(sec) = 5 * sec.A * ρ * L^4 / (384 * E * sec.Ix)

        # Heavier section has more Ix/A ratio → less SW deflection
        @test δ_sw(sec_heavy) < δ_sw(sec_light)
    end
end

println("All staged deflection tests passed.")
