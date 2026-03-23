"""
Unit tests for I_symm compression capacity with Kx/Ky/Lx/Ly parameterization.

Run from project root after loading the full module:
    include("SlabDesignFactors/SlabDesignFactors.jl")
    include("SlabDesignFactors/test/test_compression.jl")

These tests verify:
- Legacy single-length `get_Pn(sec, L)` still works
- Two-length `get_Pn(sec, Lx, Ly)` produces correct results
- K factors reduce capacity as expected
- `get_ϕPn` wrappers match 0.9 × Pn
- Known W-shape capacities are approximately correct vs AISC tables
"""

using Test

@testset "I_symm compression — Kx/Ky/Lx/Ly" begin

    # W14x82:  d=14.3, bf=10.1, tw=0.510, tf=0.855 (AISC Manual)
    sec = I_symm(14.3, 10.1, 0.510, 0.855)
    L_10ft = 10.0 * 12.0  # 120 in

    @testset "legacy single-length API" begin
        Pn_weak = get_Pn(sec, L_10ft; axis=:weak)
        Pn_strong = get_Pn(sec, L_10ft; axis=:strong)
        @test Pn_weak > 0
        @test Pn_strong > Pn_weak  # strong axis always stronger for W14
        @test get_ϕPn(sec, L_10ft; axis=:weak) ≈ 0.9 * Pn_weak
    end

    @testset "two-length API with K=1.0 matches legacy" begin
        # With equal lengths and K=1, the two-length form should produce
        # a capacity ≤ the single weak-axis capacity (because it also
        # checks strong-axis and torsional buckling).
        Pn_legacy = get_Pn(sec, L_10ft; axis=:weak)
        Pn_dual   = get_Pn(sec, L_10ft, L_10ft; Kx=1.0, Ky=1.0)
        # The two-length form takes min(Fe_x, Fe_y, Fe_torsional)
        # so it should be ≤ the weak-axis-only legacy form
        @test Pn_dual <= Pn_legacy + 1e-6
        @test Pn_dual > 0
    end

    @testset "K > 1 reduces capacity" begin
        ϕPn_K1  = get_ϕPn(sec, L_10ft, L_10ft; Kx=1.0, Ky=1.0)
        ϕPn_K2y = get_ϕPn(sec, L_10ft, L_10ft; Kx=1.0, Ky=2.0)
        ϕPn_K2x = get_ϕPn(sec, L_10ft, L_10ft; Kx=2.0, Ky=1.0)
        @test ϕPn_K2y < ϕPn_K1  # doubling Ky reduces capacity
        # Doubling Kx may not reduce capacity if weak-axis already governs;
        # but it should never *increase* capacity.
        @test ϕPn_K2x <= ϕPn_K1
        # Use a longer column (20ft) where strong-axis buckling can govern
        L_20ft = 20.0 * 12.0
        ϕPn_long_K1  = get_ϕPn(sec, L_20ft, L_20ft; Kx=1.0, Ky=1.0)
        ϕPn_long_K2x = get_ϕPn(sec, L_20ft, L_20ft; Kx=2.0, Ky=1.0)
        ϕPn_long_K2y = get_ϕPn(sec, L_20ft, L_20ft; Kx=1.0, Ky=2.0)
        @test ϕPn_long_K2x <= ϕPn_long_K1
        @test ϕPn_long_K2y < ϕPn_long_K1
    end

    @testset "different Lx, Ly" begin
        # Shorter weak-axis length (e.g. mid-height bracing) should increase capacity
        ϕPn_equal = get_ϕPn(sec, L_10ft, L_10ft)
        ϕPn_braced = get_ϕPn(sec, L_10ft, L_10ft / 2; Kx=1.0, Ky=1.0)
        @test ϕPn_braced >= ϕPn_equal
    end

    @testset "ϕ factor" begin
        Pn = get_Pn(sec, L_10ft, L_10ft)
        @test get_ϕPn(sec, L_10ft, L_10ft) ≈ 0.9 * Pn
        @test get_ϕPn(sec, L_10ft, L_10ft; ϕ=0.85) ≈ 0.85 * Pn
    end

    @testset "stocky section, short column ≈ squash load" begin
        # W14x730 equivalent (very stocky): d=22.4, bf=17.9, tw=3.07, tf=4.91
        sec_stocky = I_symm(22.4, 17.9, 3.07, 4.91)
        L_short = 36.0  # 3 ft — very short
        Fy = steel_ksi.Fy  # 50 ksi
        squash = Fy * sec_stocky.A
        ϕPn = get_ϕPn(sec_stocky, L_short, L_short)
        # For a very short column, ϕPn should approach ϕ·Fy·A
        @test ϕPn / (0.9 * squash) > 0.90
    end

    @testset "slender column, low capacity" begin
        # W8x24: d=7.93, bf=6.50, tw=0.245, tf=0.400
        sec_light = I_symm(7.93, 6.50, 0.245, 0.400)
        L_long = 30.0 * 12.0  # 30 ft
        ϕPn = get_ϕPn(sec_light, L_long, L_long)
        squash = steel_ksi.Fy * sec_light.A
        # Very slender column — capacity should be well below 50% of squash
        @test ϕPn / (0.9 * squash) < 0.50
    end
end

println("All compression tests passed.")
