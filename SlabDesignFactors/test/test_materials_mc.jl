"""
Unit tests for ConcreteMaterial scenarios and Monte Carlo machinery.

Run with:  julia SlabDesignFactors/test/test_materials_mc.jl

Self-contained — only uses base Julia + Statistics.
"""

using Test
using Statistics
using Random

# ── Bootstrap: load only the lightweight core files we need ──────────────────
include("../core/conversions.jl")
include("../core/constants.jl")
include("../core/materials.jl")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. ConcreteMaterial struct & predefined scenarios
# ═══════════════════════════════════════════════════════════════════════════════
@testset "ConcreteMaterial struct" begin
    @testset "Convenience constructor derives E_c and imperial values" begin
        mat = ConcreteMaterial("test"; ρ_concrete=2400.0, f′c_MPa=27.6, ecc_concrete=0.152)
        @test mat.name == "test"
        @test mat.ρ_concrete == 2400.0
        @test mat.f′c_MPa == 27.6
        @test mat.ecc_concrete == 0.152
        @test isapprox(mat.E_c_MPa, 4700.0 * sqrt(27.6), atol=1.0)
        @test mat.ρ_concrete_kipin3 > 0
        @test mat.f′c_ksi > 0
        @test mat.E_c_ksi > 0
        @test isapprox(mat.f′c_ksi, 27.6 / 6.89476, atol=0.01)
    end

    @testset "Four predefined scenarios exist and differ" begin
        @test length(CONCRETE_SCENARIOS) == 4
        @test haskey(CONCRETE_SCENARIOS, "NW_4ksi")
        @test haskey(CONCRETE_SCENARIOS, "NW_5ksi")
        @test haskey(CONCRETE_SCENARIOS, "LW_4ksi")
        @test haskey(CONCRETE_SCENARIOS, "LW_3ksi")

        nw4 = CONCRETE_SCENARIOS["NW_4ksi"]
        nw5 = CONCRETE_SCENARIOS["NW_5ksi"]
        lw4 = CONCRETE_SCENARIOS["LW_4ksi"]
        lw3 = CONCRETE_SCENARIOS["LW_3ksi"]

        @test nw4.ρ_concrete > lw4.ρ_concrete
        @test nw5.ρ_concrete > lw3.ρ_concrete
        @test lw4.ecc_concrete > nw4.ecc_concrete
        @test lw3.ecc_concrete > nw4.ecc_concrete
        @test nw5.E_c_ksi > nw4.E_c_ksi
        @test nw5.E_c_MPa > nw4.E_c_MPa
        @test nw4.ρ_concrete_kipin3 > lw4.ρ_concrete_kipin3
    end

    @testset "DEFAULT_CONCRETE matches NW_4ksi (legacy)" begin
        @test DEFAULT_CONCRETE.name == "NW_4ksi"
        @test DEFAULT_CONCRETE.ρ_concrete == ρ_CONCRETE
        @test DEFAULT_CONCRETE.ecc_concrete == ECC_CONCRETE
        @test isapprox(DEFAULT_CONCRETE.ρ_concrete_kipin3, ρ_CONCRETE_KIPIN3, rtol=1e-6)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Density propagation produces measurably different dead loads
# ═══════════════════════════════════════════════════════════════════════════════
@testset "Density-driven dead load differences" begin
    nw4 = CONCRETE_SCENARIOS["NW_4ksi"]
    lw4 = CONCRETE_SCENARIOS["LW_4ksi"]

    # Simulate the density-based slab dead load formula from get_scaled_model
    w_nw4 = 0.99 * nw4.ρ_concrete_kipin3 + 0.01 * ρ_STEEL_KIPIN3
    w_lw4 = 0.99 * lw4.ρ_concrete_kipin3 + 0.01 * ρ_STEEL_KIPIN3

    @test w_lw4 < w_nw4

    ratio = w_lw4 / w_nw4
    @test 0.70 < ratio < 0.85
    println("  Dead load ratio LW/NW: $(round(ratio, digits=3))  (expect ~0.77)")
end

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Monte Carlo sampling (rand_lognormal)
# ═══════════════════════════════════════════════════════════════════════════════

function rand_lognormal(rng::AbstractRNG, nominal::Real, σ_ln::Real, N::Int)
    μ_ln = log(nominal) - σ_ln^2 / 2
    return exp.(μ_ln .+ σ_ln .* randn(rng, N))
end

@testset "rand_lognormal sampling" begin
    N = 100_000
    rng = MersenneTwister(42)

    @testset "Recovers correct mean and CV" begin
        for (nominal, cv_target) in [(1.22, 0.10), (0.854, 0.10), (0.152, 0.15)]
            σ_ln = sqrt(log(1 + cv_target^2))
            samples = rand_lognormal(rng, nominal, σ_ln, N)

            sample_mean = mean(samples)
            sample_cv = std(samples) / mean(samples)

            @test isapprox(sample_mean, nominal, rtol=0.02)
            @test isapprox(sample_cv, cv_target, atol=0.015)
            println("  nominal=$nominal cv_target=$cv_target → mean=$(round(sample_mean, digits=4)) cv=$(round(sample_cv, digits=4))")
        end
    end

    @testset "Deterministic seeding" begin
        rng1 = MersenneTwister(123)
        rng2 = MersenneTwister(123)
        s1 = rand_lognormal(rng1, 1.0, 0.1, 1000)
        s2 = rand_lognormal(rng2, 1.0, 0.1, 1000)
        @test s1 == s2
    end

    @testset "All samples are positive" begin
        samples = rand_lognormal(MersenneTwister(99), 0.05, 0.3, 10_000)
        @test all(s -> s > 0, samples)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 4. mc_summary statistics
# ═══════════════════════════════════════════════════════════════════════════════

function mc_summary(samples::AbstractVector{<:Real})
    return (
        mean = mean(samples),
        std  = std(samples),
        p5   = quantile(samples, 0.05),
        p25  = quantile(samples, 0.25),
        p50  = quantile(samples, 0.50),
        p75  = quantile(samples, 0.75),
        p95  = quantile(samples, 0.95),
    )
end

@testset "mc_summary" begin
    rng = MersenneTwister(77)
    samples = randn(rng, 50_000)

    s = mc_summary(samples)
    @test isapprox(s.mean, 0.0, atol=0.02)
    @test isapprox(s.std, 1.0, atol=0.02)
    @test isapprox(s.p50, 0.0, atol=0.02)
    @test s.p5 < s.p25 < s.p50 < s.p75 < s.p95
    @test isapprox(s.p5, -1.645, atol=0.05)
    @test isapprox(s.p95, 1.645, atol=0.05)
end

# ═══════════════════════════════════════════════════════════════════════════════
# 5. End-to-end MC with synthetic structural data
# ═══════════════════════════════════════════════════════════════════════════════
@testset "End-to-end MC with synthetic masses" begin
    rng = MersenneTwister(2024)
    N = 10_000

    norm_mass_beams = 35.0   # kg/m²
    norm_mass_slab  = 280.0  # kg/m²
    norm_mass_rebar = 8.0    # kg/m²

    ecc_steel_nom = ECC_STEEL
    ecc_rebar_nom = ECC_REBAR
    ecc_conc_nom  = ECC_CONCRETE

    σ_steel = sqrt(log(1 + 0.10^2))
    σ_rebar = sqrt(log(1 + 0.10^2))
    σ_conc  = sqrt(log(1 + 0.15^2))

    ecc_steel_samples = rand_lognormal(rng, ecc_steel_nom, σ_steel, N)
    ecc_rebar_samples = rand_lognormal(rng, ecc_rebar_nom, σ_rebar, N)
    ecc_conc_samples  = rand_lognormal(rng, ecc_conc_nom,  σ_conc,  N)

    ec_beams_nom = ecc_steel_nom * norm_mass_beams
    ec_slab_nom  = ecc_conc_nom  * norm_mass_slab
    ec_rebar_nom_val = ecc_rebar_nom * norm_mass_rebar
    ec_total_nom = ec_beams_nom + ec_slab_nom + ec_rebar_nom_val

    @testset "Joint MC mean ≈ nominal total EC" begin
        ec_joint = norm_mass_beams .* ecc_steel_samples .+
                   norm_mass_slab  .* ecc_conc_samples  .+
                   norm_mass_rebar .* ecc_rebar_samples
        mc = mc_summary(ec_joint)
        @test isapprox(mc.mean, ec_total_nom, rtol=0.02)
        @test mc.std > 0
        @test mc.p5 < mc.p95
        println("  Nominal EC total: $(round(ec_total_nom, digits=2))")
        println("  MC joint mean:    $(round(mc.mean, digits=2)) ± $(round(mc.std, digits=2))")
        println("  MC joint 90% CI:  [$(round(mc.p5, digits=2)), $(round(mc.p95, digits=2))]")
    end

    @testset "Independent sweeps differ from each other" begin
        ec_steel_only = norm_mass_beams .* ecc_steel_samples .+ ec_slab_nom .+ ec_rebar_nom_val
        ec_conc_only  = ec_beams_nom .+ norm_mass_slab .* ecc_conc_samples .+ ec_rebar_nom_val
        ec_rebar_only = ec_beams_nom .+ ec_slab_nom .+ norm_mass_rebar .* ecc_rebar_samples

        s_steel = mc_summary(ec_steel_only)
        s_conc  = mc_summary(ec_conc_only)
        s_rebar = mc_summary(ec_rebar_only)

        for (label, s) in [("steel", s_steel), ("concrete", s_conc), ("rebar", s_rebar)]
            @test isapprox(s.mean, ec_total_nom, rtol=0.02)
        end

        @test s_conc.std > s_rebar.std
        println("  Steel sweep std:    $(round(s_steel.std, digits=3))")
        println("  Concrete sweep std: $(round(s_conc.std, digits=3))")
        println("  Rebar sweep std:    $(round(s_rebar.std, digits=3))")
    end

    @testset "Variance decomposition sums to ~1" begin
        v_steel = var(norm_mass_beams .* ecc_steel_samples)
        v_conc  = var(norm_mass_slab  .* ecc_conc_samples)
        v_rebar = var(norm_mass_rebar .* ecc_rebar_samples)
        v_total = v_steel + v_conc + v_rebar

        vf_steel = v_steel / v_total
        vf_conc  = v_conc  / v_total
        vf_rebar = v_rebar / v_total

        @test isapprox(vf_steel + vf_conc + vf_rebar, 1.0, atol=1e-10)
        @test vf_steel > 0
        @test vf_conc > 0
        @test vf_rebar > 0
        @test vf_conc > vf_rebar
        println("  Variance fractions: steel=$(round(vf_steel, digits=3)) conc=$(round(vf_conc, digits=3)) rebar=$(round(vf_rebar, digits=3))")
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# 6. Cross-scenario MC comparison (different concretes yield different results)
# ═══════════════════════════════════════════════════════════════════════════════
@testset "Cross-scenario MC: NW vs LW concrete" begin
    rng_nw = MersenneTwister(1000)
    rng_lw = MersenneTwister(1000)
    N = 10_000

    nw = CONCRETE_SCENARIOS["NW_4ksi"]
    lw = CONCRETE_SCENARIOS["LW_4ksi"]

    norm_mass_beams = 35.0
    norm_mass_rebar = 8.0

    slab_volume_per_m2 = 0.15  # m³/m² (150mm slab depth)
    norm_mass_slab_nw = slab_volume_per_m2 * nw.ρ_concrete * (1 - 0.01)
    norm_mass_slab_lw = slab_volume_per_m2 * lw.ρ_concrete * (1 - 0.01)

    σ_steel = sqrt(log(1 + 0.10^2))
    σ_rebar = sqrt(log(1 + 0.10^2))
    σ_conc  = sqrt(log(1 + 0.15^2))

    ec_joint_nw = norm_mass_beams .* rand_lognormal(rng_nw, ECC_STEEL, σ_steel, N) .+
                  norm_mass_slab_nw .* rand_lognormal(rng_nw, nw.ecc_concrete, σ_conc, N) .+
                  norm_mass_rebar .* rand_lognormal(rng_nw, ECC_REBAR, σ_rebar, N)

    ec_joint_lw = norm_mass_beams .* rand_lognormal(rng_lw, ECC_STEEL, σ_steel, N) .+
                  norm_mass_slab_lw .* rand_lognormal(rng_lw, lw.ecc_concrete, σ_conc, N) .+
                  norm_mass_rebar .* rand_lognormal(rng_lw, ECC_REBAR, σ_rebar, N)

    mc_nw = mc_summary(ec_joint_nw)
    mc_lw = mc_summary(ec_joint_lw)

    println("  NW_4ksi: slab_mass=$(round(norm_mass_slab_nw, digits=1)) kg/m², EC mean=$(round(mc_nw.mean, digits=2)) ± $(round(mc_nw.std, digits=2))")
    println("  LW_4ksi: slab_mass=$(round(norm_mass_slab_lw, digits=1)) kg/m², EC mean=$(round(mc_lw.mean, digits=2)) ± $(round(mc_lw.std, digits=2))")

    @test norm_mass_slab_nw > norm_mass_slab_lw
    @test nw.ecc_concrete < lw.ecc_concrete

    diff_pct = abs(mc_nw.mean - mc_lw.mean) / mc_nw.mean * 100
    @test diff_pct > 1.0
    println("  Mean EC difference: $(round(diff_pct, digits=1))%")

    @test mc_nw.p5 != mc_lw.p5
    @test mc_nw.p95 != mc_lw.p95
end

println("\n✓ All material scenario & Monte Carlo tests passed.")
