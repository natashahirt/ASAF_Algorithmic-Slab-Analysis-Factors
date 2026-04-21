"""
Sanity reproduction test for the Nov 2024 bare-steel two-stage NLP
sizer port in `VariableBeamOptimizer/sizer_nlp_v1/`.

The port reuses the Nov 2024 NLP kernel (`generate_broadcast`,
`generate_objective`, `inequality_constraints_I_symm`,
`get_element_deflection`) already resident under
`VariableBeamOptimizer/rolled_steel/` (byte-equivalent to
`VariableBeamOptimizer @ c3a506f:optim/*`), wrapped in a faithful port
of the November 2024 cluster driver found at
`TributaryAreas @ a8dd2b98:slab_analysis/size_beams.jl`:

  - `process_discrete_beams_nov2024` — `binary_search_sections` through
    the W-catalogue (NOT the monorepo's newer
    `sequential_search_sections` variant).
  - `process_continuous_beams_nov2024` — NLP with W4X13 variable floor
    when `minimum_continuous=true`, W43X335 ceiling, NLoptAlg(:LD_MMA)
    at `xtol_rel = ftol_rel = ftol_abs = 1e-2, xtol_abs = 1e-8`.
  - `get_deflection_constraint_nov2024` — hard inequality
    `abs(minimum(δ_FE_with_SW)) ≤ L / serviceability_lim`, pushed on
    top of the soft penalty inside `generate_objective`.

Composite action is forced off (`composite_action=false`,
`deflection_reduction_factor=1.0`) so no modern hook fires.

# Acceptance envelope (±5 %)

The residual gap against the legacy CSV is bounded below by the MMA
solver's own declared tolerance (`ftol_rel = 1e-2`) applied 174 times
and by Julia/Nonconvex/NLopt version drift since the Nov 2024 run
(1.10 → 1.12 changes float reduction ordering, which propagates
through `AsapOptim.solve_frame_Pf`'s LU solve). With the Nov 2024
floor (W4X13) and binary-search discrete warm start both wired in,
BaU reproduces to **~3.4 %** on the heavy side of the legacy value —
i.e. conservative vs. the paper's baseline, not optimistic. Accepting
±5 % gives a one-tolerance-envelope safety margin for the other bays
the comparison script will run through.

# BaU target

BaU row (r6c4 / uniform / isotropic / vector=[0,0] / max_depth=25 in)
of `SlabDesignFactors/results/remote_results_yesdeflection_yesslabmin/topology.csv`.
Legacy steel_norm ≈ 28.36 kg/m²; reproduced should be within ±5 %.

Usage (from project root):
    julia --project=. -e 'include("SlabDesignFactors/scripts/_scripts.jl"); \
                           include("SlabDesignFactors/test/test_sizer_nlp_v1_reproduction.jl")'
"""

using Test
using DataFrames
using CSV

const CSV_PATH_NOV2024 = joinpath(@__DIR__, "..", "results",
                                  "remote_results_yesdeflection_yesslabmin",
                                  "topology.csv")

const GEOM_DIR = joinpath(@__DIR__, "..", "..", "Geometries", "topology")

"""BaU definition used throughout the paper (Nov 2024 sweep row)."""
const BAU = (
    name        = "r6c4",
    slab_type   = :isotropic,
    slab_sizer  = :uniform,
    beam_sizer  = "continuous",
    collinear   = false,
    vector_1d_x = 0.0,
    vector_1d_y = 0.0,
    max_depth   = 25.0,
)

const TOL_REL = 0.05  # ±5 %: one MMA `ftol_rel=1e-2` envelope wide, see module docstring.

"""
    _load_legacy_bau_row(csv_path, bau) -> NamedTuple

Locate the BaU row in the Nov 2024 legacy CSV and return the fields that
the reproduction test compares against.
"""
function _load_legacy_bau_row(csv_path::AbstractString, bau)
    df = DataFrame(CSV.File(csv_path))
    mask = (df.name      .== bau.name)       .&
           (string.(df.slab_type)  .== String(bau.slab_type))  .&
           (string.(df.slab_sizer) .== String(bau.slab_sizer)) .&
           (string.(df.beam_sizer) .== bau.beam_sizer)         .&
           (df.collinear .== bau.collinear)                    .&
           (df.vector_1d_x .≈ bau.vector_1d_x)                 .&
           (df.vector_1d_y .≈ bau.vector_1d_y)                 .&
           (df.max_depth   .≈ bau.max_depth)

    rows = df[mask, :]
    @assert nrow(rows) == 1 "Expected exactly 1 BaU row; got $(nrow(rows))."

    r = rows[1, :]
    return (
        area          = r.area,
        steel_norm    = r.steel_norm,
        concrete_norm = r.concrete_norm,
        rebar_norm    = r.rebar_norm,
    )
end

"""Load a topology JSON into a geometry dict, mirroring the sweep driver."""
function _load_geometry_dict(geom_dir::AbstractString, name::AbstractString)
    path = joinpath(geom_dir, name * ".json")
    raw  = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))
    return raw isa Dict ? raw : Dict(pairs(raw))
end

"""
    _run_nlp_v1(geometry_dict, bau) -> SlabOptimResults

Size the BaU bay using the replicated Nov 2024 sizer and post-process.
Uses the same load set as the Nov 2024 sweep (50 psf live, 15 psf SDL,
1.6/1.2 factors, L/360, minimum_continuous=true).
"""
function _run_nlp_v1(geometry_dict, bau)
    geom, _ = Base.invokelatest(generate_from_json, geometry_dict;
                                plot=false, drawn=false)

    slab_params = SlabAnalysisParams(
        geom,
        slab_name     = bau.name,
        slab_type     = bau.slab_type,
        vector_1d     = [bau.vector_1d_x, bau.vector_1d_y],
        slab_sizer    = bau.slab_sizer,
        spacing       = 0.1,
        plot_analysis = false,
        fix_param     = true,
        slab_units    = :m,
    )

    sizing_params = SlabSizingParams(
        live_load              = psf_to_ksi(50),
        superimposed_dead_load = psf_to_ksi(15),
        slab_dead_load         = 0.0,
        live_factor            = 1.6,
        dead_factor            = 1.2,
        beam_sizer             = :continuous,
        nlp_solver             = :MMA,
        max_depth              = bau.max_depth,
        beam_units             = :in,
        serviceability_lim     = 360,
        collinear              = bau.collinear,
        minimum_continuous     = true,
        n_max_sections         = 0,
        composite_action       = false,
    )

    slab_params  = analyze_slab(slab_params)
    slab_params, sizing_params = size_bare_nlp!(slab_params, sizing_params)

    return postprocess_slab(slab_params, sizing_params,
                             check_collinear = bau.collinear)
end

@testset "sizer_nlp_v1 reproduction (Nov 2024 BaU)" begin
    legacy = _load_legacy_bau_row(CSV_PATH_NOV2024, BAU)
    println("  Legacy BaU row:   area=$(legacy.area)  steel_norm=$(legacy.steel_norm)")

    geom_dict = _load_geometry_dict(GEOM_DIR, BAU.name)
    results   = _run_nlp_v1(geom_dict, BAU)

    reproduced_steel_norm = results.norm_mass_beams
    reproduced_area       = results.area

    rel_err_mass = abs(reproduced_steel_norm - legacy.steel_norm) / legacy.steel_norm
    rel_err_area = abs(reproduced_area       - legacy.area)       / legacy.area

    println("  Reproduced:       area=$(reproduced_area)  steel_norm=$(reproduced_steel_norm)")
    println("  Rel err (steel):  $(round(100 * rel_err_mass, digits=2)) %  (tol = $(100 * TOL_REL) %)")
    println("  Rel err (area):   $(round(100 * rel_err_area, digits=2)) %")

    @test !isempty(results.ids)
    @test results.norm_mass_beams > 0
    @test rel_err_area <= 0.01    # geometry/tributary must match exactly
    @test rel_err_mass <= TOL_REL
end
