"""
    process_continuous_beams_nov2024(params::SlabSizingParams; initial_vars=[]) -> SlabSizingParams

Replicates the Nov 2024 continuous-NLP beam sizer as it existed at
VariableBeamOptimizer @ c3a506f (`optim/optimize_I_symm.jl` +
`optim/utils.jl` + `setup/problem/beam_problem.jl`), for validation against
today's composite-aware pipeline. The function body mirrors
`sizer/size_beams.jl::process_continuous_beams` with two deliberate
constraints so the Nov 2024 regime is reproduced exactly:

  1. Composite action is forced off (`composite_action = false`,
     `deflection_reduction_factor = 1.0`) regardless of caller settings.
     Nov 2024 had no composite-aware branch.
  2. The hard deflection constraint is the **bare-steel, non-composite
     fallback** branch of `get_deflection_constraint` — i.e. the
     `composite_action=false` path that enforces
     `(δ_point + δ_beam_sw) ≤ L / serviceability_lim` with
     `deflection_reduction_factor = 1.0`. This matches the Nov 2024
     regime: strict bare-steel L/360 on the bare section, serialised by
     `get_element_deflection(..., dead_load=0.0)` plus an analytical
     self-weight term. The soft `(δ - δ_max)^2 / (δ - δ_max)^4` penalty
     inside `generate_objective` remains in play; it is numerically
     dominated by the hard constraint near the feasible boundary.

Every other numerical choice — the Nonconvex/NLoptAlg(:LD_MMA) solver,
`xtol_rel=1e-2, ftol_abs=1e-2, ftol_rel=1e-2` tolerances, the
`[0.01, 0.01, 0.001, 0.001]` vs `W6X8.5` minimum-variable bounds, and the
`W43X335` maximum bounds — is preserved from `c3a506f`.

This routine is sequential (no demand iteration, no staged loop). It is
intended for bay-by-bay validation, not production sweeps.

# Arguments
- `params::SlabSizingParams`: Sizing parameters. Must have a model with
    `params.model.elements[:beam]` populated and `params.load_dictionary`
    constructed for the same model (typically done by `optimal_beamsizer`
    preamble; `size_bare_nlp!` handles this automatically).
- `initial_vars`: Optional warm-start sections (same conventions as the
    legacy driver: `Vector{String}` of `W_imperial` names, a
    `Vector{Vector{Float64}}` of `[h, w, tw, tf]` per beam, or a
    `Vector{Section}`). Defaults to `W6X8.5` per beam.

# Returns
The input `params`, mutated with `minimizers`, `minimums`, `ids`,
`M_maxs`, `V_maxs`, `x_maxs` populated.
"""
function process_continuous_beams_nov2024(params::SlabSizingParams; initial_vars::Vector=Any[])

    params.composite_action            = false
    params.deflection_reduction_factor = 1.0

    model         = params.model
    beam_elements = model.elements[:beam]

    for i in 1:lastindex(beam_elements)

        println("=== $i / $(lastindex(beam_elements)) === (nov2024 NLP)")

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

        optimization_model = Nonconvex.Model()

        fixed_variables = [true, true, true, true]
        broadcast_in, broadcast_out = generate_broadcast(fixed_variables, n_sections)

        if params.minimum_continuous
            min_h, min_w, min_tw, min_tf = get_geometry_vars(W_imperial("W4X13"))
        else
            min_h, min_w, min_tw, min_tf = [0.01, 0.01, 0.001, 0.001]
        end

        max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))

        if !iszero(params.max_beam_depth) && !isinf(params.max_beam_depth)
            max_h = min(params.max_beam_depth, max_h)
        end

        addvar!(optimization_model,
                broadcast_in([min_h, min_w, min_tw, min_tf]),
                broadcast_in([max_h, max_w, max_tw, max_tf]))

        beam_params = get_frameoptparams(params, beam_elements[i],
                                         params.load_dictionary[get_element_id(beam_elements[i])])

        objective = generate_objective(params, beam_params, beam_forces,
                                       broadcast_out, i_check, n_sections,
                                       interpolation=akima_interpolation)
        set_objective!(optimization_model, objective)

        constraints = inequality_constraints_I_symm(beam_forces, n_sections,
                                                    i_check=i_check,
                                                    interpolation_function=akima_interpolation,
                                                    broadcast_out=broadcast_out)

        if params.deflection_limit
            push!(constraints,
                  get_deflection_constraint_nov2024(beam_params, beam_forces.x[end], params))
        end

        for constraint in constraints
            add_ineq_constraint!(optimization_model, constraint)
        end

        alg = NLoptAlg(:LD_MMA)
        options = NLoptOptions(
            xtol_rel = 1e-2,
            xtol_abs = 1e-8,
            ftol_abs = 1e-2,
            ftol_rel = 1e-2,
        )

        if isempty(initial_vars)
            init_vars = broadcast_in(get_geometry_vars(W_imperial("W4X13")))
        elseif eltype(initial_vars) <: AbstractString
            init_vars = broadcast_in(get_geometry_vars(W_imperial(initial_vars[i])))
        elseif eltype(initial_vars) <: AbstractVector
            init_vars = broadcast_in(Vector{Float64}(initial_vars[i]))
        elseif eltype(initial_vars) <: Asap.AbstractSection
            init_vars = broadcast_in(get_geometry_vars(initial_vars[i]))
        else
            init_vars = broadcast_in(get_geometry_vars(W_imperial("W4X13")))
        end

        result = Nonconvex.optimize(optimization_model, alg, init_vars, options=options)

        println("")
        push!(params.minimizers, result.minimizer)
        push!(params.minimums,   result.minimum)
        push!(params.ids,        string(round(I_symm(result.minimizer...).A, digits=2)))
        push!(params.M_maxs,     M_max)
        push!(params.V_maxs,     V_max)
        push!(params.x_maxs,     x_max)
    end

    return params
end
