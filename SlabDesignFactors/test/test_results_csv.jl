"""
CSV export contract:
- hard failures (`result_ok == false`) must not populate CSV with masses, sections,
  or utilizations that could be read as a valid design.
- serviceability warnings remain exportable so downstream filtering can exclude them later.
"""

using Test

@testset "create_results_dataframe omits design quantities when result_ok is false" begin
    r = SlabOptimResults(
        slab_name         = "fail_case",
        slab_type         = :isotropic,
        vector_1d         = [1.0, 0.0],
        slab_sizer        = :uniform,
        beam_sizer        = :discrete,
        area              = 120.0,
        max_depth         = 40.0,
        norm_mass_beams   = 99.0,
        norm_mass_columns = 88.0,
        norm_mass_slab    = 77.0,
        norm_mass_rebar   = 66.0,
        norm_mass_fireproofing = 5.0,
        sections          = ["W21X44", "W18X35"],
        ids               = ["b1", "b2"],
        nlp_solver        = "MIP",
        deflection_limit  = true,
        collinear         = true,
        max_bay_span      = 360.0,
        max_util_M        = 1.5,
        max_util_V        = 0.9,
        max_col_util      = 0.5,
        geometry_file     = "fail_case",
        span_ok           = true,
        result_ok         = false,
        strength_ok       = false,
        serviceability_ok = true,
        staged_ok         = true,
        column_ok         = true,
        solver_status     = "STRENGTH_FAIL",
        diagnostic_flags  = "strength_fail;result_not_ok",
        diagnostic_messages = "Triggered flags: strength_fail, result_not_ok",
        config_hash       = "abc123",
    )
    df = create_results_dataframe([r], false)
    @test df.steel_norm[1] == 0.0
    @test df.column_norm[1] == 0.0
    @test df.concrete_norm[1] == 0.0
    @test df.rebar_norm[1] == 0.0
    @test df.fireproofing_norm[1] == 0.0
    @test df.sections[1] == "Any[]"
    @test df.ids[1] == "Any[]"
    @test df.max_util_M[1] == 0.0
    @test df.max_util_V[1] == 0.0
    @test df.result_ok[1] == false
    @test df.span_ok[1] == true
    @test df.strength_ok[1] == false
    @test df.staged_deflection_limit[1] == true
    @test df.solver_status[1] == "STRENGTH_FAIL"
    @test occursin("strength_fail", df.diagnostic_flags[1])
    @test df.config_hash[1] == "abc123"
end

@testset "create_results_dataframe preserves design quantities for serviceability warnings" begin
    r = SlabOptimResults(
        slab_name         = "warn_case",
        slab_type         = :isotropic,
        vector_1d         = [1.0, 0.0],
        slab_sizer        = :uniform,
        beam_sizer        = :discrete,
        area              = 120.0,
        max_depth         = 40.0,
        norm_mass_beams   = 99.0,
        norm_mass_columns = 8.0,
        norm_mass_slab    = 77.0,
        norm_mass_rebar   = 6.0,
        norm_mass_fireproofing = 5.0,
        sections          = ["W21X44", "W18X35"],
        ids               = ["b1", "b2"],
        col_sections      = ["W10X12"],
        col_Pu            = [10.0],
        col_ϕPn           = [20.0],
        col_util          = [0.5],
        δ_slab_dead       = [0.1, 0.2],
        δ_beam_dead       = [0.01, 0.02],
        δ_sdl             = [0.03, 0.04],
        δ_live            = [0.05, 0.06],
        δ_total           = [0.19, 0.32],
        Δ_limit_live      = [0.2, 0.2],
        Δ_limit_total     = [0.25, 0.25],
        δ_live_ok         = [true, true],
        δ_total_ok        = [true, false],
        composite_action  = true,
        staged_converged  = false,
        staged_n_violations = 1,
        staged_ok         = false,
        n_L360_fail       = 0,
        n_L240_fail       = 1,
        i_L360_fail       = Int[],
        i_L240_fail       = [2],
        max_δ_total       = 0.32,
        max_bay_span      = 360.0,
        global_δ_ok       = true,
        max_util_M        = 0.95,
        max_util_V        = 0.6,
        max_col_util      = 0.5,
        geometry_file     = "warn_case",
        span_ok           = true,
        result_ok         = true,
        strength_ok       = true,
        serviceability_ok = false,
        column_ok         = true,
        solver_status     = "SERVICEABILITY_WARNING",
        diagnostic_flags  = "serviceability_fail;l240_fail;staged_not_converged;staged_violations_remaining",
        diagnostic_messages = "Triggered flags: serviceability_fail, l240_fail, staged_not_converged, staged_violations_remaining",
        nlp_solver        = "MIP",
        deflection_limit  = true,
        collinear         = true,
        config_hash       = "warn123",
    )
    df = create_results_dataframe([r], false)
    @test df.steel_norm[1] == 99.0
    @test df.sections[1] != "Any[]"
    @test df.result_ok[1] == true
    @test df.span_ok[1] == true
    @test df.strength_ok[1] == true
    @test df.serviceability_ok[1] == false
    @test df.staged_deflection_limit[1] == true
    @test df.staged_ok[1] == false
    @test df.solver_status[1] == "SERVICEABILITY_WARNING"
    @test occursin("serviceability_fail", df.diagnostic_flags[1])
end
