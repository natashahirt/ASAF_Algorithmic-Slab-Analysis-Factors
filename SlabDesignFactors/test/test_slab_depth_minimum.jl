"""
Unit tests for `SlabAnalysisParams.slab_depth_minimum` and span-based depth clamping.

Requires project load via `SlabDesignFactors/test/run.jl` (same as other unit tests).
"""

using Test

const _TOPO_JSON = joinpath(@__DIR__, "..", "..", "Geometries", "topology", "r1c1.json")

function _r1c1_geom()
    d = geometry_dict_from_json_path(_TOPO_JSON)
    geom, _ = generate_from_json(d; plot=false, drawn=false)
    return geom
end

@testset "slab_depth_minimum constructor" begin
    geom = _r1c1_geom()
    p = SlabAnalysisParams(geom; slab_depth_minimum=1e-9, slab_units=:m)
    @test p.slab_depth_minimum ≈ 0.001
    p2 = SlabAnalysisParams(geom; slab_depth_minimum=0.0, slab_units=:m)
    @test p2.slab_depth_minimum ≈ 0.001
    p_default = SlabAnalysisParams(geom; slab_units=:m)
    @test p_default.slab_depth_minimum ≈ 0.001
    p3 = SlabAnalysisParams(geom; slab_depth_minimum=0.125, slab_units=:m)
    @test p3.slab_depth_minimum ≈ 0.125
    # 0.001 m = 1 mm when working in mm
    pmm = SlabAnalysisParams(geom; slab_depth_minimum=0.0005, slab_units=:mm)
    @test pmm.slab_depth_minimum ≈ 1.0
    @test_throws AssertionError SlabAnalysisParams(geom; slab_depth_minimum=-0.01, slab_units=:m)
end

@testset "slab_depth_minimum clamp in get_slab_depths (via analyze_slab)" begin
    geom = _r1c1_geom()
    base = SlabAnalysisParams(
        geom;
        slab_name="t",
        slab_type=:isotropic,
        vector_1d=[1.0, 0.0],
        slab_sizer=:uniform,
        spacing=0.1,
        plot_analysis=false,
        fix_param=true,
        slab_units=:m,
        slab_depth_minimum=0.0,
    )
    hi = SlabAnalysisParams(
        geom;
        slab_name="t",
        slab_type=:isotropic,
        vector_1d=[1.0, 0.0],
        slab_sizer=:uniform,
        spacing=0.1,
        plot_analysis=false,
        fix_param=true,
        slab_units=:m,
        slab_depth_minimum=0.5,
    )
    base = analyze_slab(base)
    hi = analyze_slab(hi)
    @test !isempty(hi.slab_depths)
    @test all(d >= 0.5 - 1e-9 for d in hi.slab_depths)
    @test all(hi.slab_depths .>= base.slab_depths .- 1e-9)
end
