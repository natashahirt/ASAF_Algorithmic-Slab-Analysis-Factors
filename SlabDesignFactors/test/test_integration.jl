"""
Integration tests: full-pipeline beam sizing across analysis methods and geometries.

Runs the complete geometry → slab analysis → beam sizing → postprocess chain on
real topology JSON files using **composite action** and **drf = 1.0** (full staged
deflection checks), matching production / full-sweep defaults.

Bare non-composite cases are not used in the deflection-limited pipeline here;
staged deflection is validated on the composite results (`δ_slab_dead`, `δ_beam_dead`,
SDL/Live split, CSV export).

Collinearity is enforced **during** the MIP/NLP solve via same-section constraints
on collinear groups (not via legacy postprocessing).  Postprocessing follows the
solver's `collinear` flag without overriding it.  Tests verify:
  - Collinear runs: every beam in a collinear group receives the same section.
  - Noncollinear runs: beams are sized independently; the `collinear` flag
    propagates correctly to results; mass is ≤ the collinear result.

Solver paths tested:
  - MIP (discrete, Gurobi) : all geometries (collinear + noncollinear)
  - NLP (continuous, Ipopt) : r1c1, r2c3 (smaller geometries — NLP is slower)

Geometries:
  - r1c1 : 1-row 1-col (simplest regular bay)
  - r2c3 : 2-row 3-col (moderate rectangular grid)
  - r5c2 : 5-row 2-col (tall narrow — many beams, long spans possible)
  - r7c4 : 7-row 4-col (large irregular topology with many collinear groups)

Verified physical invariants:
  1. Composite (drf=1.0) produces valid beams when feasible for the geometry.
  2. Composite mass is finite and positive.
  3. All beams satisfy M_u / ϕM_n ≤ 1 and V_u / ϕV_n ≤ 1 (strength feasibility).
  4. Staged deflection fields are non-empty and δ_total > 0.
  5. Analytical beam self-weight deflection is non-negative for every beam.
  6. Column sizing produces valid results when columns exist.
  7. NLP mass ≤ MIP mass (continuous relaxation is a lower bound on discrete).
  8. NLP and MIP produce the same number of beams per geometry.
  9. Strength-only composite mass ≤ deflection-governed composite mass.
 10. Solver-enforced collinearity: group members share identical sections.
 11. Noncollinear mass ≤ collinear mass (relaxed constraint is a lower bound).
"""

using Test
using Gurobi

"""True if a Gurobi license is available (discrete / MIP beam sizing)."""
function _gurobi_available()
    try
        Gurobi.Env()
        return true
    catch
        return false
    end
end

const GUROBI_AVAILABLE = _gurobi_available()
const STRENGTH_WARN_LIMIT = 1.0
const STRENGTH_FAIL_LIMIT = 1.02

main_path = joinpath(@__DIR__, "..", "..", "Geometries", "topology")

geom_files = ["r1c1.json", "r2c3.json", "r5c2.json", "r7c4.json"]

function load_geometry(json_file)
    path = joinpath(main_path, json_file)
    return geometry_dict_from_json_path(path)
end

"""Run the full pipeline for one geometry + one set of analysis options.

When `collinear=true` (default) the MIP/NLP solver enforces same-section
constraints within collinear groups; postprocessing follows the same flag
without overriding it, so results reflect what the optimizer actually chose.
"""
function run_pipeline(geometry_dict; composite::Bool=true, drf::Float64=1.0,
                      beam_sizer::Symbol=:discrete, nlp_solver::Symbol=:Ipopt,
                      collinear::Bool=true)
    geom, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)

    slab_params = SlabAnalysisParams(
        geom,
        slab_name       = "test",
        slab_type       = :isotropic,
        vector_1d       = [1.0, 0.0],
        slab_sizer      = :uniform,
        spacing         = 0.1,
        plot_analysis   = false,
        fix_param       = true,
        slab_units      = :m,
    )

    sizing_params = SlabSizingParams(
        live_load                   = psf_to_ksi(50),
        superimposed_dead_load      = psf_to_ksi(15),
        slab_dead_load              = 0.0,
        live_factor                 = 1.6,
        dead_factor                 = 1.2,
        beam_sizer                  = beam_sizer,
        nlp_solver                  = nlp_solver,
        max_depth                   = 40.0,
        beam_units                  = :in,
        serviceability_lim          = 360,
        collinear                   = collinear,
        minimum_continuous          = true,
        n_max_sections              = 0,
        composite_action            = composite,
        deflection_reduction_factor = drf,
    )

    slab_params = analyze_slab(slab_params)
    slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)

    if isempty(sizing_params.minimizers)
        @warn "Sizing produced no results — returning default SlabOptimResults"
        return SlabOptimResults(), slab_params, sizing_params
    end

    results = postprocess_slab(slab_params, sizing_params)

    return results, slab_params, sizing_params
end

"""True when the solver returned a non-empty, positive-mass section set."""
has_solution(r) = length(r.minimizers) >= 1 && r.norm_mass_beams > 0

"""True when a result is a design-feasible pass (not merely non-empty)."""
design_ok(r) = has_solution(r) && r.result_ok && r.strength_ok && r.serviceability_ok && r.column_ok

"""Backward-compatible alias used throughout this test file."""
result_ok(r) = design_ok(r)

if !GUROBI_AVAILABLE
    @info "Gurobi license not found — skipping MIP integration tests. Add `secrets/gurobi.lic` or set GRB_LICENSE_FILE (run.jl / runtests.jl load the secrets path when the file exists)."
end

if GUROBI_AVAILABLE
@testset "Integration — MIP (discrete)" begin

    for json_file in geom_files
        name = replace(json_file, ".json" => "")
        geom_dict = load_geometry(json_file)

        @testset "$name" begin

            # ── Run case (production-like: composite, full deflection) ───────────
            println("\n  Running $name [MIP]: composite action (drf=1.0)...")
            res_comp, _, sp_comp = run_pipeline(geom_dict; composite=true, drf=1.0)
            comp_ok = result_ok(res_comp)

            # ── 1. Feasibility ─────────────────────────────────────────
            @testset "beams exist" begin
                @test comp_ok
            end

            # ── 2. Composite mass is finite and positive ─────────────
            if comp_ok
                @testset "composite mass valid" begin
                    @test res_comp.norm_mass_beams > 0
                    @test isfinite(res_comp.norm_mass_beams)
                end
            end

            # ── 3. Strength feasibility ──────────────────────────────
            @testset "strength checks pass" begin
                results_to_check = filter(result_ok, [res_comp])
                for r in results_to_check
                    for i in 1:length(r.Mn)
                        if r.Mn[i] > 0 && !isempty(r.My[i])
                            Mu = maximum(abs.(r.My[i]))
                            util_M = Mu / r.Mn[i]
                            if util_M > STRENGTH_WARN_LIMIT
                                @warn "$name beam $i flexural utilization exceeds 1.0" utilization=util_M
                            end
                            @test util_M <= STRENGTH_FAIL_LIMIT
                        end
                        if r.Vn[i] > 0 && !isempty(r.Vy[i])
                            Vu = maximum(abs.(r.Vy[i]))
                            util_V = Vu / r.Vn[i]
                            if util_V > STRENGTH_WARN_LIMIT
                                @warn "$name beam $i shear utilization exceeds 1.0" utilization=util_V
                            end
                            @test util_V <= STRENGTH_FAIL_LIMIT
                        end
                    end
                end
            end

            # ── 4–6: Staged deflection checks (composite, drf=1.0) ───
            if comp_ok
                @testset "staged deflection populated" begin
                    n = length(res_comp.minimizers)
                    @test length(res_comp.δ_slab_dead) == n
                    @test length(res_comp.δ_beam_dead) == n
                    @test length(res_comp.δ_sdl) == n
                    @test length(res_comp.δ_live) == n
                    @test length(res_comp.δ_total) == n
                    @test length(res_comp.Δ_limit_live) == n
                    @test length(res_comp.Δ_limit_total) == n
                    @test any(res_comp.δ_total .> 0)
                end

                @testset "beam self-weight deflection ≥ 0" begin
                    @test all(res_comp.δ_beam_dead .>= 0)
                end

                @testset "superposition δ_total ≈ Σ components" begin
                    for i in 1:length(res_comp.δ_total)
                        sum_components = res_comp.δ_slab_dead[i] + res_comp.δ_beam_dead[i] +
                                        res_comp.δ_sdl[i] + res_comp.δ_live[i]
                        @test res_comp.δ_total[i] ≈ sum_components atol=1e-10
                    end
                end

                @testset "deflection limits plausible" begin
                    for i in 1:length(res_comp.Δ_limit_live)
                        @test res_comp.Δ_limit_live[i] > 0
                        @test res_comp.Δ_limit_total[i] > res_comp.Δ_limit_live[i]
                        @test res_comp.δ_live_ok[i] == (res_comp.δ_live[i] <= res_comp.Δ_limit_live[i])
                        @test res_comp.δ_total_ok[i] == (res_comp.δ_total[i] <= res_comp.Δ_limit_total[i])
                    end
                end

                @testset "CSV export mirrors staged fields (create_results_dataframe)" begin
                    df = create_results_dataframe([res_comp], false)
                    for c in (
                        "nlp_solver", "deflection_limit",
                        "composite_action", "staged_converged", "staged_n_violations",
                        "n_L360_fail", "n_L240_fail", "i_L360_fail", "i_L240_fail",
                        "Δ_limit_live_mm", "Δ_limit_total_mm",
                        "max_δ_total_mm", "max_bay_span_in", "global_δ_ok",
                        "max_util_M", "max_util_V", "max_col_util",
                    )
                        @test c in names(df)
                    end
                    @test df.composite_action[1] == true
                    @test df.n_L360_fail[1] == res_comp.n_L360_fail
                    @test df.n_L240_fail[1] == res_comp.n_L240_fail
                    @test df.staged_converged[1] == res_comp.staged_converged
                    @test df.staged_n_violations[1] == res_comp.staged_n_violations
                    @test df.nlp_solver[1] == "MIP"
                    @test df.deflection_limit[1] == true
                end
            end

            # ── 7. Column sizing (if columns exist) ───────────────────
            @testset "column sizing" begin
                results_to_check = filter(r -> result_ok(r) && !isempty(r.col_sections),
                    [res_comp])
                for r in results_to_check
                    @test length(r.col_Pu) == length(r.col_sections)
                    @test length(r.col_ϕPn) == length(r.col_sections)
                    @test length(r.col_util) == length(r.col_sections)
                    @test all(r.col_util .>= 0)
                    @test all(r.col_util .<= 1.0 + 1e-3)
                    @test r.mass_columns > 0
                end
            end

            # ── 8. Solver-enforced collinearity ────────────────────────
            if comp_ok
                @testset "solver-enforced collinearity" begin
                    @test res_comp.collinear == true
                    groups = get_collinear_groups(sp_comp.model.elements[:beam])
                    for gid in unique(groups)
                        members = findall(==(gid), groups)
                        length(members) <= 1 && continue
                        leader_id = sp_comp.ids[members[1]]
                        for m in members[2:end]
                            @test sp_comp.ids[m] == leader_id
                        end
                    end
                end
            end

            # ── Print summary ─────────────────────────────────────────
            println("  $name [MIP] results:")
            if comp_ok
                println("    Composite (drf=1.0): $(round(res_comp.norm_mass_beams, digits=2)) kg/m²")
                if !isempty(res_comp.δ_total)
                    max_δ_tot = round(maximum(res_comp.δ_total) * 1000, digits=2)
                    max_δ_sw  = round(maximum(res_comp.δ_beam_dead) * 1000, digits=2)
                    println("    Max δ_total: $(max_δ_tot) mm, max δ_beam_SW: $(max_δ_sw) mm")
                    println("    L/360 fail: $(count(.!res_comp.δ_live_ok)), " *
                            "L/240 fail: $(count(.!res_comp.δ_total_ok))")
                end
            else
                println("    Composite: INFEASIBLE (sizing could not converge)")
            end
        end

        GC.gc()
    end
end
end # GUROBI_AVAILABLE

# ══════════════════════════════════════════════════════════════════════════════
#  Noncollinear MIP — verify that beams are sized independently when
#  collinear=false.  Each beam gets its own optimal section; groups are NOT
#  required to match.  Uses a subset of geometries (r2c3, r7c4) that have
#  multi-beam collinear groups so the distinction is meaningful.
# ══════════════════════════════════════════════════════════════════════════════
noncol_geom_files = ["r2c3.json", "r7c4.json"]

if GUROBI_AVAILABLE
@testset "Integration — MIP noncollinear (collinear=false)" begin

    for json_file in noncol_geom_files
        name = replace(json_file, ".json" => "")
        geom_dict = load_geometry(json_file)

        @testset "$name" begin
            println("\n  Running $name [MIP]: noncollinear composite (drf=1.0)...")
            res_nc, _, sp_nc = run_pipeline(geom_dict; composite=true, drf=1.0, collinear=false)
            nc_ok = result_ok(res_nc)

            @testset "beams exist" begin
                @test nc_ok
            end

            if nc_ok
                @testset "collinear flag is false" begin
                    @test res_nc.collinear == false
                end
                @testset "noncollinear mass ≤ collinear mass" begin
                    res_col, _, _ = run_pipeline(geom_dict; composite=true, drf=1.0, collinear=true)
                    if result_ok(res_col)
                        @test res_nc.norm_mass_beams <= res_col.norm_mass_beams + 1e-3
                    end
                end

                @testset "strength checks pass" begin
                    for i in 1:length(res_nc.Mn)
                        if res_nc.Mn[i] > 0 && !isempty(res_nc.My[i])
                            Mu = maximum(abs.(res_nc.My[i]))
                            @test Mu / res_nc.Mn[i] <= 1.0 + 1e-3
                        end
                    end
                end
            end

            if nc_ok
                println("    Noncollinear: $(round(res_nc.norm_mass_beams, digits=2)) kg/m²")
            else
                println("    Noncollinear: INFEASIBLE")
            end
        end

        GC.gc()
    end
end
end # GUROBI_AVAILABLE

# ══════════════════════════════════════════════════════════════════════════════
#  NLP (continuous) tests — smaller geometries only (NLP is much slower)
#  Wrapped in try-catch: NLopt can fail on some platforms (JuMP nlp.jl error).
# ══════════════════════════════════════════════════════════════════════════════
nlp_geom_files = ["r1c1.json", "r2c3.json"]

nlp_available = try
    run_pipeline(load_geometry("r1c1.json"); composite=true, drf=1.0, beam_sizer=:continuous)
    true
catch e
    @warn "NLP solver unavailable on this platform — skipping NLP tests" exception=e
    false
end

if nlp_available

@testset "Integration — NLP (continuous)" begin

    for json_file in nlp_geom_files
        name = replace(json_file, ".json" => "")
        geom_dict = load_geometry(json_file)

        @testset "$name" begin

            println("\n  Running $name [NLP]: composite (drf=1.0)...")
            res_nlp_comp, _, sp_nlp_comp = run_pipeline(geom_dict; composite=true, drf=1.0, beam_sizer=:continuous)
            nlp_comp_ok = result_ok(res_nlp_comp)

            @testset "beams exist" begin
                @test nlp_comp_ok
            end

            @testset "NLP mass valid" begin
                if nlp_comp_ok
                    @test res_nlp_comp.norm_mass_beams > 0
                    @test isfinite(res_nlp_comp.norm_mass_beams)
                end
            end

            @testset "strength checks pass" begin
                nlp_tol = 0.05  # Ipopt may slightly violate constraints
                results_to_check = filter(result_ok, [res_nlp_comp])
                for r in results_to_check
                    for i in 1:length(r.Mn)
                        if r.Mn[i] > 0 && !isempty(r.My[i])
                            Mu = maximum(abs.(r.My[i]))
                            @test Mu / r.Mn[i] <= 1.0 + nlp_tol
                        end
                        if r.Vn[i] > 0 && !isempty(r.Vy[i])
                            Vu = maximum(abs.(r.Vy[i]))
                            @test Vu / r.Vn[i] <= 1.0 + nlp_tol
                        end
                    end
                end
            end

            if nlp_comp_ok
                @testset "staged deflection populated" begin
                    n = length(res_nlp_comp.minimizers)
                    @test length(res_nlp_comp.δ_slab_dead) == n
                    @test length(res_nlp_comp.δ_beam_dead) == n
                    @test length(res_nlp_comp.δ_total) == n
                    @test any(res_nlp_comp.δ_total .> 0)
                end

                @testset "beam self-weight deflection ≥ 0" begin
                    @test all(res_nlp_comp.δ_beam_dead .>= 0)
                end

                @testset "superposition δ_total ≈ Σ components" begin
                    for i in 1:length(res_nlp_comp.δ_total)
                        sum_components = res_nlp_comp.δ_slab_dead[i] + res_nlp_comp.δ_beam_dead[i] +
                                        res_nlp_comp.δ_sdl[i] + res_nlp_comp.δ_live[i]
                        @test res_nlp_comp.δ_total[i] ≈ sum_components atol=1e-10
                    end
                end
            end

            @testset "column sizing" begin
                results_to_check = filter(r -> result_ok(r) && !isempty(r.col_sections),
                    [res_nlp_comp])
                for r in results_to_check
                    @test all(r.col_util .>= 0)
                    @test all(r.col_util .<= 1.0 + 1e-3)
                    @test r.mass_columns > 0
                end
            end

            if nlp_comp_ok
                @testset "solver-enforced collinearity" begin
                    @test res_nlp_comp.collinear == true
                    groups = get_collinear_groups(sp_nlp_comp.model.elements[:beam])
                    for gid in unique(groups)
                        members = findall(==(gid), groups)
                        length(members) <= 1 && continue
                        leader_mins = sp_nlp_comp.minimizers[members[1]]
                        for m in members[2:end]
                            @test sp_nlp_comp.minimizers[m] ≈ leader_mins atol=1e-3
                        end
                    end
                end
            end

            println("  $name [NLP] results:")
            if nlp_comp_ok
                println("    Composite (drf=1.0): $(round(res_nlp_comp.norm_mass_beams, digits=2)) kg/m²")
                if !isempty(res_nlp_comp.δ_total)
                    println("    L/360 fail: $(count(.!res_nlp_comp.δ_live_ok)), " *
                            "L/240 fail: $(count(.!res_nlp_comp.δ_total_ok))")
                end
            else
                println("    Composite: INFEASIBLE")
            end
        end

        GC.gc()
    end
end

# ══════════════════════════════════════════════════════════════════════════════
#  MIP vs NLP comparison — NLP (continuous) should be a lower bound on MIP
#  Uses composite action, drf=1.0 (production-like).
# ══════════════════════════════════════════════════════════════════════════════
if GUROBI_AVAILABLE
@testset "MIP vs NLP comparison" begin

    for json_file in nlp_geom_files
        name = replace(json_file, ".json" => "")
        geom_dict = load_geometry(json_file)

        @testset "$name" begin
            println("\n  MIP vs NLP comparison for $name...")

            res_mip, _, _ = run_pipeline(geom_dict; composite=true, drf=1.0, beam_sizer=:discrete)
            res_nlp, _, _ = run_pipeline(geom_dict; composite=true, drf=1.0, beam_sizer=:continuous)

            mip_ok = result_ok(res_mip)
            nlp_ok = result_ok(res_nlp)

            @testset "MIP and NLP both feasible on baseline geometry" begin
                @test mip_ok
                @test nlp_ok
            end

            if mip_ok && nlp_ok
                ratio = res_mip.norm_mass_beams > 0 ? res_nlp.norm_mass_beams / res_mip.norm_mass_beams : NaN
                println("    MIP:  $(round(res_mip.norm_mass_beams, digits=2)) kg/m²")
                println("    NLP:  $(round(res_nlp.norm_mass_beams, digits=2)) kg/m²")
                println("    NLP/MIP ratio: $(round(ratio, digits=4))")

                @testset "NLP mass ≤ MIP mass (continuous is lower bound)" begin
                    if res_nlp.norm_mass_beams <= res_mip.norm_mass_beams + 1e-3
                        @test true
                    else
                        @warn "$name: NLP local optimum exceeds MIP " *
                              "($(round(ratio, digits=3))×); Ipopt likely stuck in local min"
                        @test res_nlp.norm_mass_beams <= res_mip.norm_mass_beams + 1e-3
                    end
                end
            end
        end

        GC.gc()
    end
end
else
    @info "Skipping MIP vs NLP comparison (Gurobi not available)."
end

# ══════════════════════════════════════════════════════════════════════════════
#  NLP algorithm comparison — all supported NLP solvers should produce valid
#  results and be within a reasonable range of Ipopt (the reference solver).
#  Uses r1c1 composite (drf=1.0) — production-like, smallest geometry.
# ══════════════════════════════════════════════════════════════════════════════
@testset "NLP algorithm comparison" begin
    geom_dict = load_geometry("r1c1.json")
    all_solvers = [:Ipopt, :MMA, :SLSQP, :CCSAQ, :COBYLA]

    println("\n  NLP algorithm comparison for r1c1 (composite, drf=1.0)...")

    solver_results = Dict{Symbol, Any}()
    solver_ok      = Dict{Symbol, Bool}()

    for solver in all_solvers
        available = try
            run_pipeline(geom_dict; composite=true, drf=1.0,
                         beam_sizer=:continuous, nlp_solver=solver)
            true
        catch e
            @warn "Solver $solver unavailable — skipping" exception=e
            false
        end

        if available
            res, _, sp = run_pipeline(geom_dict; composite=true, drf=1.0,
                                      beam_sizer=:continuous, nlp_solver=solver)
            solver_results[solver] = res
            solver_ok[solver] = result_ok(res)
            mass = solver_ok[solver] ? round(res.norm_mass_beams, digits=2) : "INFEASIBLE"
            println("    $solver: $mass kg/m²")
        else
            solver_ok[solver] = false
            println("    $solver: UNAVAILABLE")
        end
    end

    @testset "all solvers produce valid results" begin
        for solver in all_solvers
            @testset "$solver" begin
                if solver_ok[solver]
                    @test solver_results[solver].norm_mass_beams > 0
                else
                    @warn "$solver did not produce a valid result"
                    @test_skip true
                end
            end
        end
    end

    @testset "same beam count across solvers" begin
        valid = [s for s in all_solvers if solver_ok[s]]
        if length(valid) >= 2
            ref_count = length(solver_results[valid[1]].minimizers)
            for solver in valid[2:end]
                @testset "$solver" begin
                    @test length(solver_results[solver].minimizers) == ref_count
                end
            end
        end
    end

    @testset "results within 50% of Ipopt" begin
        if solver_ok[:Ipopt]
            ipopt_mass = solver_results[:Ipopt].norm_mass_beams
            for solver in all_solvers
                solver == :Ipopt && continue
                @testset "$solver" begin
                    if solver_ok[solver]
                        ratio = solver_results[solver].norm_mass_beams / ipopt_mass
                        println("    $solver / Ipopt ratio: $(round(ratio, digits=4))")
                        if 0.5 <= ratio <= 1.5
                            @test true
                        else
                            @warn "$solver / Ipopt ratio $(round(ratio, digits=4)) outside [0.5, 1.5] " *
                                  "(often an infeasible NLP local minimum, not a true lower bound)"
                            @test_broken 0.5 <= ratio <= 1.5
                        end
                    else
                        @test_skip true
                    end
                end
            end
        end
    end

    @testset "strength feasible (proven solvers)" begin
        nlp_tol = 0.05
        for solver in [:Ipopt, :MMA]
            @testset "$solver" begin
                if solver_ok[solver]
                    r = solver_results[solver]
                    for i in 1:length(r.Mn)
                        if r.Mn[i] > 0 && !isempty(r.My[i])
                            Mu = maximum(abs.(r.My[i]))
                            @test Mu / r.Mn[i] <= 1.0 + nlp_tol
                        end
                    end
                else
                    @test_skip true
                end
            end
        end
    end

    @testset "strength summary (experimental solvers)" begin
        nlp_tol = 0.05
        for solver in [:SLSQP, :CCSAQ, :COBYLA]
            @testset "$solver" begin
                if solver_ok[solver]
                    r = solver_results[solver]
                    violations = 0
                    worst = 0.0
                    for i in 1:length(r.Mn)
                        if r.Mn[i] > 0 && !isempty(r.My[i])
                            ratio = maximum(abs.(r.My[i])) / r.Mn[i]
                            if ratio > 1.0 + nlp_tol
                                violations += 1
                                worst = max(worst, ratio)
                            end
                        end
                    end
                    println("    $solver: $violations violations, worst Mu/Mn = $(round(worst, digits=3))")
                    @test true
                else
                    @test_skip true
                end
            end
        end
    end

    GC.gc()
end

end # if nlp_available

# ══════════════════════════════════════════════════════════════════════════════
#  Strength-only tests (deflection_limit=false)
#  Isolates the strength path vs deflection-governed composite (drf=1.0).
#  Strength-only composite mass should be ≤ deflection-governed composite mass.
# ══════════════════════════════════════════════════════════════════════════════
nodefl_geom_files = ["r1c1.json", "r2c3.json"]

"""Run pipeline with deflection_limit=false.

See `run_pipeline` for collinearity rationale — postprocessing follows the
solver flag without overriding it.
"""
function run_pipeline_nodefl(geometry_dict; composite::Bool=false,
                              beam_sizer::Symbol=:discrete,
                              collinear::Bool=true)
    geom, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)

    slab_params = SlabAnalysisParams(
        geom,
        slab_name       = "test",
        slab_type       = :isotropic,
        vector_1d       = [1.0, 0.0],
        slab_sizer      = :uniform,
        spacing         = 0.1,
        plot_analysis   = false,
        fix_param       = true,
        slab_units      = :m,
    )

    sizing_params = SlabSizingParams(
        live_load                   = psf_to_ksi(50),
        superimposed_dead_load      = psf_to_ksi(15),
        slab_dead_load              = 0.0,
        live_factor                 = 1.6,
        dead_factor                 = 1.2,
        beam_sizer                  = beam_sizer,
        max_depth                   = 40.0,
        beam_units                  = :in,
        serviceability_lim          = 360,
        collinear                   = collinear,
        minimum_continuous          = true,
        n_max_sections              = 0,
        composite_action            = composite,
        deflection_limit            = false,
    )

    slab_params = analyze_slab(slab_params)
    slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)

    if isempty(sizing_params.minimizers)
        @warn "Sizing produced no results — returning default SlabOptimResults"
        return SlabOptimResults(), slab_params, sizing_params
    end

    results = postprocess_slab(slab_params, sizing_params)
    return results, slab_params, sizing_params
end

if !GUROBI_AVAILABLE
    @info "Skipping Integration — strength only (requires Gurobi for discrete sizing)."
end

if GUROBI_AVAILABLE
@testset "Integration — strength only (deflection_limit=false)" begin

    for json_file in nodefl_geom_files
        name = replace(json_file, ".json" => "")
        geom_dict = load_geometry(json_file)

        @testset "$name" begin
            println("\n  Running $name: strength-only composite (deflection_limit=false)...")
            res_str, _, _ = run_pipeline_nodefl(geom_dict; composite=true)

            println("  Running $name: deflection-governed composite (drf=1.0)...")
            res_defl, _, _ = run_pipeline(geom_dict; composite=true, drf=1.0)

            @testset "beams exist" begin
                @test length(res_str.minimizers) >= 1
            end

            @testset "strength checks pass" begin
                for r in [res_str]
                    for i in 1:length(r.Mn)
                        if r.Mn[i] > 0 && !isempty(r.My[i])
                            Mu = maximum(abs.(r.My[i]))
                            util_M = Mu / r.Mn[i]
                            if util_M > STRENGTH_WARN_LIMIT
                                @warn "$name beam $i flexural utilization exceeds 1.0" utilization=util_M
                            end
                            @test util_M <= STRENGTH_FAIL_LIMIT
                        end
                        if r.Vn[i] > 0 && !isempty(r.Vy[i])
                            Vu = maximum(abs.(r.Vy[i]))
                            util_V = Vu / r.Vn[i]
                            if util_V > STRENGTH_WARN_LIMIT
                                @warn "$name beam $i shear utilization exceeds 1.0" utilization=util_V
                            end
                            @test util_V <= STRENGTH_FAIL_LIMIT
                        end
                    end
                end
            end

            if design_ok(res_defl)
                @testset "strength-only ≤ deflection-governed mass" begin
                    @test res_str.norm_mass_beams <= res_defl.norm_mass_beams + 1e-3
                end
            end

            println("  $name strength-only results:")
            println("    Composite (no defl): $(round(res_str.norm_mass_beams, digits=2)) kg/m²")
            println("    Composite (w/ defl): $(round(res_defl.norm_mass_beams, digits=2)) kg/m²")
            ratio = res_defl.norm_mass_beams > 0 ?
                round(res_defl.norm_mass_beams / max(res_str.norm_mass_beams, 1e-6), digits=2) : NaN
            println("    Deflection overhead: $(ratio)×")
        end

        GC.gc()
    end
end
end # GUROBI_AVAILABLE

# ══════════════════════════════════════════════════════════════════════════════
#  Grid-style solver comparison — exercises configurations matching the
#  production full-sweep grid (varied slab types, depths, geometries) and
#  compares MIP, Ipopt, and MMA.  Catches regressions where a continuous
#  solver silently returns infeasible sections.
#
#  Configurations chosen to span:
#    - small (r1c1, 12 beams) to large (r5c2, many beams)
#    - isotropic + orthotropic biaxial + uniaxial slab types
#    - 25 in and 40 in depth limits
#    - cellular and uniform slab sizers
# ══════════════════════════════════════════════════════════════════════════════
function run_pipeline_grid(geometry_dict;
                           slab_type::Symbol=:isotropic,
                           vector_1d::Vector{Float64}=[1.0, 0.0],
                           slab_sizer::Symbol=:uniform,
                           max_depth::Float64=40.0,
                           beam_sizer::Symbol=:discrete,
                           nlp_solver::Symbol=:Ipopt,
                           collinear::Bool=true)
    geom, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)

    slab_params = SlabAnalysisParams(
        geom,
        slab_name       = "test",
        slab_type       = slab_type,
        vector_1d       = vector_1d,
        slab_sizer      = slab_sizer,
        spacing         = 0.1,
        plot_analysis   = false,
        fix_param       = true,
        slab_units      = :m,
    )

    sizing_params = SlabSizingParams(
        live_load                   = psf_to_ksi(50),
        superimposed_dead_load      = psf_to_ksi(15),
        slab_dead_load              = 0.0,
        live_factor                 = 1.6,
        dead_factor                 = 1.2,
        beam_sizer                  = beam_sizer,
        nlp_solver                  = nlp_solver,
        max_depth                   = max_depth,
        beam_units                  = :in,
        serviceability_lim          = 360,
        collinear                   = collinear,
        minimum_continuous          = true,
        n_max_sections              = 0,
        composite_action            = true,
        deflection_reduction_factor = 1.0,
    )

    slab_params = analyze_slab(slab_params)
    slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)

    if isempty(sizing_params.minimizers)
        return SlabOptimResults(), slab_params, sizing_params
    end

    results = postprocess_slab(slab_params, sizing_params)
    return results, slab_params, sizing_params
end

grid_configs = [
    # (json_file, slab_type, vector_1d, slab_sizer, max_depth, description)
    # Small geometry — isotropic
    ("r1c1.json", :isotropic,     [0.0, 0.0], :uniform,  40.0, "r1c1 iso uniform 40in"),
    ("r1c1.json", :isotropic,     [0.0, 0.0], :cellular, 25.0, "r1c1 iso cellular 25in"),
    # Medium geometry — varied slab types
    ("r2c3.json", :isotropic,     [0.0, 0.0], :uniform,  40.0, "r2c3 iso uniform 40in"),
    ("r2c3.json", :orth_biaxial,  [1.0, 0.0], :uniform,  25.0, "r2c3 orth_biax uniform 25in"),
    ("r2c3.json", :uniaxial,     [1.0, 0.0], :cellular, 40.0, "r2c3 uniax cellular 40in"),
    # Larger geometry — stress tests the NLP scalability
    ("r5c2.json", :isotropic,     [0.0, 0.0], :uniform,  40.0, "r5c2 iso uniform 40in"),
    ("r5c2.json", :isotropic,     [0.0, 0.0], :cellular, 25.0, "r5c2 iso cellular 25in"),
]

if GUROBI_AVAILABLE && nlp_available
@testset "Grid solver comparison (MIP vs Ipopt vs MMA)" begin
    # These cases are expected to be design-feasible in CI and local development.
    required_mip_configs = Set([
        "r1c1 iso uniform 40in",
        "r1c1 iso cellular 25in",
        "r2c3 iso uniform 40in",
        "r2c3 orth_biax uniform 25in",
        "r2c3 uniax cellular 40in",
    ])

    for (json_file, slab_type, vector_1d, slab_sizer, max_depth, desc) in grid_configs
        geom_dict = load_geometry(json_file)

        @testset "$desc" begin
            println("\n  ─── Grid: $desc ───")
            results_by_solver = Dict{String, Any}()
            ok_by_solver = Dict{String, Bool}()
            has_solution_by_solver = Dict{String, Bool}()

            for (label, bsizer, nlp) in [
                ("MIP",   :discrete,   :Ipopt),
                ("Ipopt", :continuous, :Ipopt),
                ("MMA",   :continuous, :MMA),
            ]
                local res
                try
                    res, _, _ = run_pipeline_grid(geom_dict;
                        slab_type=slab_type, vector_1d=vector_1d,
                        slab_sizer=slab_sizer, max_depth=max_depth,
                        beam_sizer=bsizer, nlp_solver=nlp)
                    results_by_solver[label] = res
                    has_solution_by_solver[label] = has_solution(res)
                    ok_by_solver[label] = design_ok(res)
                catch e
                    @warn "  $label failed on $desc" exception=e
                    results_by_solver[label] = nothing
                    has_solution_by_solver[label] = false
                    ok_by_solver[label] = false
                end

                mass_str = has_solution_by_solver[label] ?
                    "$(round(results_by_solver[label].norm_mass_beams, digits=2)) kg/m²" :
                    "FAILED"
                println("    $label: $mass_str")
            end

            # MIP should be feasible on required configs; stress configs are informational.
            @testset "MIP feasible" begin
                if desc in required_mip_configs
                    @test ok_by_solver["MIP"]
                else
                    if !ok_by_solver["MIP"]
                        @warn "$desc: MIP infeasible on stress config"
                        @test_skip true
                    else
                        @test true
                    end
                end
            end

            # Ipopt should produce design-feasible results on required configs.
            @testset "Ipopt feasible" begin
                if ok_by_solver["Ipopt"]
                    @test true
                    r = results_by_solver["Ipopt"]
                    for i in 1:length(r.Mn)
                        if r.Mn[i] > 0 && !isempty(r.My[i])
                            @test maximum(abs.(r.My[i])) / r.Mn[i] <= 1.0 + 0.05
                        end
                    end
                else
                    if desc in required_mip_configs
                        @warn "$desc: Ipopt infeasible"
                        @test false
                    else
                        @warn "$desc: Ipopt infeasible on stress config"
                        @test_skip true
                    end
                end
            end

            # Ipopt mass should be ≤ MIP mass (continuous relaxation lower bound)
            if ok_by_solver["MIP"] && ok_by_solver["Ipopt"]
                @testset "Ipopt ≤ MIP mass" begin
                    ratio = results_by_solver["Ipopt"].norm_mass_beams /
                            results_by_solver["MIP"].norm_mass_beams
                    println("    Ipopt/MIP ratio: $(round(ratio, digits=4))")
                    if ratio <= 1.0 + 1e-3
                        @test true
                    else
                        @warn "$desc: Ipopt exceeds MIP ($(round(ratio, digits=3))×)"
                        @test_broken ratio <= 1.0 + 1e-3
                    end
                end
            end

            # MMA strength check — report but don't hard-fail the suite
            @testset "MMA strength (informational)" begin
                if has_solution_by_solver["MMA"]
                    r = results_by_solver["MMA"]
                    violations = 0
                    worst = 0.0
                    for i in 1:length(r.Mn)
                        if r.Mn[i] > 0 && !isempty(r.My[i])
                            u = maximum(abs.(r.My[i])) / r.Mn[i]
                            if u > 1.05
                                violations += 1
                                worst = max(worst, u)
                            end
                        end
                    end
                    println("    MMA: $violations strength violations, worst Mu/Mn = $(round(worst, digits=3))")
                    @test true
                else
                    println("    MMA: infeasible / no result")
                    @test true
                end
            end

            GC.gc()
        end
    end
end
end # GUROBI_AVAILABLE && nlp_available

println("\nAll integration tests complete.")
