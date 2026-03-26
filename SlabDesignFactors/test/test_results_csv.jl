"""
CSV export contract: failed designs (`result_ok == false`) must not populate CSV with
masses, sections, or utilizations that could be read as a valid design.
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
        result_ok         = false,
        strength_ok       = false,
        serviceability_ok = true,
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
    @test df.strength_ok[1] == false
    @test df.solver_status[1] == "STRENGTH_FAIL"
    @test occursin("strength_fail", df.diagnostic_flags[1])
    @test df.config_hash[1] == "abc123"
end
