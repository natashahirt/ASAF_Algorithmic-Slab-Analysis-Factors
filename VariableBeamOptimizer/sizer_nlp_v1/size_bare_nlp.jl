"""
    size_bare_nlp!(slab_params::SlabAnalysisParams, sizing_params::SlabSizingParams;
                  initial_vars=[], serviceability_lim=360) -> (SlabAnalysisParams, SlabSizingParams)

Reproduce the Nov 2024 bare-steel two-stage sizing workflow exactly:
discrete W-catalog sequential search first, then continuous NLP
warm-started from the discrete sections.

This mirrors `VariableBeamOptimizer @ c43446d::sizer/size_beams.jl::iterate_discrete_continuous`,
which is the earliest public ancestor of the Nov 24 2024 sweep that
produced `SlabDesignFactors/results/remote_results_*/topology.csv`
(the sweep itself ran with VBO pinned at c3a506f plus an unpublished
driver that was first released in c43446d). The two-stage nature is
essential: the continuous NLP uses MMA with loose `ftol=1e-2`
tolerances and is sensitive to the warm start; seeding from the
discrete section converges to the correct global basin, whereas
cold-starting from `W6X8.5` leaves MMA at a much lighter local
optimum.

Settings pinned here regardless of what the caller set, so the Nov
2024 regime is reproduced bit-for-bit:

  - `beam_sizer`                  = `:continuous` (after discrete pass)
  - `composite_action`            = `false` (Nov 2024 had no composite
                                              branch; forces today's
                                              `get_deflection_constraint`
                                              down its bare-steel
                                              fallback, which matches
                                              c43446d's simpler
                                              3-arg formula once
                                              `deflection_reduction_factor = 1.0`).
  - `deflection_limit`            = `true`
  - `deflection_reduction_factor` = `1.0`
  - `serviceability_lim`          = `360` (overrideable for edge cases).

Post-processing, collinear grouping, and CSV emission are left to
the caller via `postprocess_slab`; this facade only fills in
`sizing_params.minimizers`, `.ids`, `.M_maxs`, `.V_maxs`,
`.x_maxs`, and `.minimums`.
"""
function size_bare_nlp!(slab_params::SlabAnalysisParams, sizing_params::SlabSizingParams;
                        initial_vars::Vector=Any[],
                        serviceability_lim::Real=360)

    sizing_params.composite_action            = false
    sizing_params.deflection_limit            = true
    sizing_params.deflection_reduction_factor = 1.0
    sizing_params.serviceability_lim          = serviceability_lim

    sizing_params.model = slab_params.model

    if !isempty(sizing_params.M_maxs)
        slab_params   = reset_SlabAnalysisParams(slab_params, slab_params.model)
        sizing_params = reset_SlabSizingParams(sizing_params)
    end

    conversion_factor = convert_to_m[slab_params.slab_units] * 1 / convert_to_m[sizing_params.beam_units]
    sizing_params.area = slab_params.area * conversion_factor^2

    if sizing_params.max_assembly_depth
        slab_depth                   = maximum(slab_params.slab_depths) * conversion_factor
        sizing_params.max_beam_depth = sizing_params.max_depth - slab_depth
    end

    sizing_params.max_bay_span = maximum(slab_params.max_spans) * conversion_factor

    sizing_params.model           = get_scaled_model(slab_params, sizing_params, conversion_factor)
    sizing_params.load_dictionary = get_load_dictionary_by_id(sizing_params.model)

    if sizing_params.collinear === true
        sizing_params.collinear_groups = get_collinear_groups(sizing_params.model.elements[:beam])
    end

    # ── Stage 1: discrete binary-search warm start ────────────────────────
    # TributaryAreas @ a8dd2b98 used `binary_search_sections` (not the
    # monorepo's current `sequential_search_sections`) for the discrete
    # stage. `process_discrete_beams_nov2024` ports that verbatim.
    # Harvest its minimizers as the NLP warm start.
    sizing_params.beam_sizer = :discrete
    try
        sizing_params = process_discrete_beams_nov2024(sizing_params)
    catch e
        if e isa NoValidSectionsError
            @warn "Nov 2024 discrete warm-start produced no feasible sections" exception=e
            sizing_params.minimizers = Vector{Float64}[]
        else
            rethrow(e)
        end
    end

    warm_start_vars = copy(sizing_params.minimizers)

    # ── Stage 2: continuous NLP warm-started from discrete ────────────────
    sizing_params.beam_sizer = :continuous
    sizing_params            = reset_SlabSizingParams(sizing_params)
    sizing_params.composite_action            = false
    sizing_params.deflection_limit            = true
    sizing_params.deflection_reduction_factor = 1.0
    sizing_params.serviceability_lim          = serviceability_lim

    continuous_init = isempty(initial_vars) ? warm_start_vars : initial_vars

    try
        sizing_params = process_continuous_beams_nov2024(sizing_params; initial_vars=continuous_init)
    catch e
        if e isa NoValidSectionsError
            @warn "Nov 2024 NLP sizing produced no feasible sections" exception=e
            sizing_params.minimizers = Vector{Float64}[]
        else
            rethrow(e)
        end
    end

    sizing_params.staged_converged    = true
    sizing_params.staged_n_violations = 0

    return slab_params, sizing_params
end
