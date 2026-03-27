"""
Fast preflight parity checks for baseline geometries.

Purpose:
- Fail fast before long integration/full-suite runs.
- Ensure MIP and NLP both produce design-feasible results on baseline cases.
- Catch staged-deflection and strength regressions in the MIP/NLP handoff.
"""

using Test
using Gurobi

function _gurobi_available_preflight()
    try
        Gurobi.Env()
        return true
    catch
        return false
    end
end

const _PREFLIGHT_GUROBI = _gurobi_available_preflight()
const _PREFLIGHT_GEOMS = ["r1c1.json", "r2c3.json"]
const _PREFLIGHT_PATH = joinpath(@__DIR__, "..", "..", "Geometries", "topology")
const _PREFLIGHT_STRENGTH_LIMIT = 1.02

function _preflight_run_pipeline(json_file::String; beam_sizer::Symbol)
    geom_dict = geometry_dict_from_json_path(joinpath(_PREFLIGHT_PATH, json_file))
    geom, _ = Base.invokelatest(generate_from_json, geom_dict; plot=false, drawn=false)

    slab_params = SlabAnalysisParams(
        geom,
        slab_name="preflight",
        slab_type=:isotropic,
        vector_1d=[1.0, 0.0],
        slab_sizer=:uniform,
        spacing=0.1,
        plot_analysis=false,
        fix_param=true,
        slab_units=:m,
    )

    sizing_params = SlabSizingParams(
        live_load=psf_to_ksi(50),
        superimposed_dead_load=psf_to_ksi(15),
        slab_dead_load=0.0,
        live_factor=1.6,
        dead_factor=1.2,
        beam_sizer=beam_sizer,
        nlp_solver=:Ipopt,
        max_depth=40.0,
        beam_units=:in,
        serviceability_lim=360,
        collinear=true,
        minimum_continuous=true,
        n_max_sections=0,
        composite_action=true,
        deflection_reduction_factor=1.0,
    )

    slab_params = analyze_slab(slab_params)
    slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)

    if isempty(sizing_params.minimizers)
        return SlabOptimResults(), slab_params, sizing_params
    end
    results = postprocess_slab(slab_params, sizing_params)
    return results, slab_params, sizing_params
end

_design_ok_preflight(r) = !isempty(r.minimizers) && r.result_ok && r.strength_ok && r.serviceability_ok && r.column_ok

if _PREFLIGHT_GUROBI
@testset "Preflight MIP/NLP baseline parity" begin
    for json_file in _PREFLIGHT_GEOMS
        name = replace(json_file, ".json" => "")
        @testset "$name" begin
            res_mip, _, sp_mip = _preflight_run_pipeline(json_file; beam_sizer=:discrete)
            res_nlp, _, sp_nlp = _preflight_run_pipeline(json_file; beam_sizer=:continuous)

            @test _design_ok_preflight(res_mip)
            @test _design_ok_preflight(res_nlp)
            @test res_mip.staged_n_violations == 0
            @test res_nlp.staged_n_violations == 0

            @test res_mip.max_util_M <= _PREFLIGHT_STRENGTH_LIMIT
            @test res_mip.max_util_V <= _PREFLIGHT_STRENGTH_LIMIT
            @test res_nlp.max_util_M <= _PREFLIGHT_STRENGTH_LIMIT
            @test res_nlp.max_util_V <= _PREFLIGHT_STRENGTH_LIMIT

            @test length(res_mip.minimizers) == length(res_nlp.minimizers)
            @test res_nlp.norm_mass_beams <= res_mip.norm_mass_beams + 1e-3

            @testset "stored final demands match postprocess" begin
                for (res, sp) in ((res_mip, sp_mip), (res_nlp, sp_nlp))
                    @test length(sp.M_maxs) == length(res.My)
                    for i in 1:length(sp.M_maxs)
                        if !isempty(res.My[i])
                            @test isapprox(
                                sp.M_maxs[i],
                                maximum(abs.(res.My[i]));
                                rtol=5e-2, atol=1e-3,
                            )
                        end
                    end
                end
            end
        end
    end
end
else
    @info "Skipping preflight MIP/NLP parity (Gurobi unavailable)."
end

