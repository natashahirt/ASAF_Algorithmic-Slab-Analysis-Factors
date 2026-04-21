"""
    process_discrete_beams_nov2024(params::SlabSizingParams) -> SlabSizingParams

Replicates the Nov 2024 discrete W-catalog sizer as it existed at
`TributaryAreas @ a8dd2b98:slab_analysis/size_beams.jl::process_discrete_beams`
(the authoritative source that fed the Nov 24 2024 sweep; the later
`VariableBeamOptimizer @ c43446d` driver was a near-verbatim re-home
of this same routine).

Differences from the monorepo's current `process_discrete_beams`:

  1. Uses `binary_search_sections` (logarithmic probe of a
     weight-sorted W-catalogue, returns the lightest feasible shape)
     rather than `sequential_search_sections` (linear scan). On a
     sorted catalogue the two agree; on the unsorted catalogue
     Nov 2024 actually used, they can pick different ties and that
     shifts the warm start the NLP sees.
  2. Calls `get_deflection_constraint` with the simple 3-arg form so
     today's `composite_action=false` branch is driven without the
     analytical self-weight correction path ever firing — matching the
     Nov 2024 regime exactly.

Returns the mutated `params` with `minimizers`, `ids`, `minimums`,
`M_maxs`, `V_maxs`, `x_maxs` populated. Beams that hit `check_uniqueness!`
are still filled via the cached prior result.
"""
function process_discrete_beams_nov2024(params::SlabSizingParams)

    model         = params.model
    beam_elements = model.elements[:beam]

    for i in 1:lastindex(beam_elements)

        println("=== $i / $(lastindex(beam_elements)) === (nov2024 discrete)")

        n_sections  = 1
        i_check     = Int64[]
        beam_loads  = params.load_dictionary[get_element_id(beam_elements[i])]
        beam_forces = InternalForces(beam_elements[i], beam_loads, resolution=200)

        M_max = maximum(abs.(beam_forces.My))
        V_max = maximum(abs.(beam_forces.Vy))
        x_max = maximum(abs.(beam_forces.x))

        params, unique = check_uniqueness!(params, M_max, V_max, x_max)
        if !unique
            continue
        end

        fixed_variables = [true, true, true, true]
        broadcast_in, broadcast_out = generate_broadcast(fixed_variables, n_sections)

        max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))

        if !iszero(params.max_beam_depth) && !isinf(params.max_beam_depth)
            max_h = min(params.max_beam_depth, max_h)
        end

        beam_params = get_frameoptparams(params, beam_elements[i],
                                         params.load_dictionary[get_element_id(beam_elements[i])])

        objective   = generate_objective(params, beam_params, beam_forces,
                                         broadcast_out, i_check, n_sections,
                                         interpolation=akima_interpolation)
        constraints = inequality_constraints_I_symm(beam_forces, n_sections,
                                                    i_check=i_check,
                                                    interpolation_function=akima_interpolation,
                                                    broadcast_out=broadcast_out)

        if params.deflection_limit
            push!(constraints,
                  get_deflection_constraint_nov2024(beam_params, beam_forces.x[end], params))
        end

        minimum_val, minimum_section = binary_search_sections(params.catalog_discrete,
                                                              objective, constraints, max_h)

        minimizer = get_geometry_vars(minimum_section)
        id        = minimum_section.name

        println("")
        push!(params.minimizers, minimizer)
        push!(params.minimums,   minimum_val)
        push!(params.ids,        id)
        push!(params.M_maxs,     M_max)
        push!(params.V_maxs,     V_max)
        push!(params.x_maxs,     x_max)
    end

    return params
end
