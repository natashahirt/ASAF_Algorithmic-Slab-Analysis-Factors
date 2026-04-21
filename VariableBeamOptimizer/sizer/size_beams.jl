"""
    optimal_beamsizer!(self::SlabAnalysisParams, initial_vars::Vector; max_depth::Real=21, sizing_unit::Symbol=:in, deflection_limit::Bool=true, verbose::Bool=true, minimum::Bool=true, max_assembly_depth::Bool=true)

Optimizes beam sizes based on initial variables and parameters. Returns the optimal beam sizes.
"""
function optimal_beamsizer!(self::SlabAnalysisParams, params::SlabSizingParams; initial_vars::Vector=[])
    @assert !isempty(initial_vars) "You need to input initial variables as a list of strings ['W6X8.5', ...] or vectors of floats [[4.17, 4.055, 0.27, 0.34], ...]"
    return optimal_beamsizer(self, params, initial_vars=initial_vars)
end

"""
    optimal_beamsizer(self, params; initial_vars)

Size beams under a single-pass, single-limit regime.

# Regime

Both composite (`params.composite_action = true`) and bare-steel
(`false`) modes use a **single deflection constraint** on the total
unfactored service load (see `get_deflection_constraint`). Composite
mode evaluates that constraint against the AISC I3.1a transformed-section
`I_x` (no construction-stage bookkeeping, no separate L/360-live and
L/240-total checks). Bare-steel mode uses bare `I_x` with
`params.deflection_reduction_factor` for back-compatibility with
Nov 2024-style designs that approximated composite action via a
stiffness multiplier rather than a transformed section.

No outer convergence loop, no Gurobi MIP warm start, no staged demand
iteration. If you want to warm-start the continuous NLP from a discrete
catalogue pass, run this function once with `beam_sizer = :discrete`
and feed its minimizers back as `initial_vars` to a second call with
`beam_sizer = :continuous` (see the `iterate_discrete_continuous`
pattern below).

# Arguments

- `self`   : `SlabAnalysisParams` with loads, geometry, model already
            populated by `analyze_slab`.
- `params` : `SlabSizingParams` controlling the sizing regime.
- `initial_vars`: warm-start sections for the continuous NLP (empty
            vector, `Vector{String}` of W-imperial names, a
            `Vector{Vector{Float64}}` of `[h, w, tw, tf]` per beam, or a
            `Vector{Section}`).

# Returns

`(self, params)` with `params.minimizers`, `params.minimums`,
`params.ids`, `params.M_maxs`, `params.V_maxs`, `params.x_maxs`
populated.
"""
function optimal_beamsizer(self::SlabAnalysisParams, params::SlabSizingParams;
                           initial_vars::Vector=[])

    params.model = self.model

    if !isempty(params.M_maxs)
        self   = reset_SlabAnalysisParams(self, self.model)
        params = reset_SlabSizingParams(params)
    end

    conversion_factor = convert_to_m[self.slab_units] * 1 / convert_to_m[params.beam_units]
    params.area       = self.area * conversion_factor^2

    if params.max_assembly_depth
        slab_depth            = maximum(self.slab_depths) * conversion_factor
        params.max_beam_depth = params.max_depth - slab_depth
    end

    # Composite-specific context: slab depth and perimeter set drive the
    # transformed-Ix calculation in `get_deflection_constraint`.
    if params.composite_action
        params.slab_depth_in = maximum(self.slab_depths) * conversion_factor
        params.i_perimeter   = Set(self.i_perimeter)
    end
    params.max_bay_span = maximum(self.max_spans) * conversion_factor

    params.model           = get_scaled_model(self, params, conversion_factor)
    params.load_dictionary = get_load_dictionary_by_id(params.model)

    if params.collinear == true
        params.collinear_groups = get_collinear_groups(params.model.elements[:beam])
    end

    # Serial sizers are routed through `get_deflection_constraint` directly
    # — the parallel/MIP paths in `size_fast.jl` reimplement staged deflection
    # logic internally and are incompatible with the single-limit regime.
    try
        if params.beam_sizer == :continuous
            params = process_continuous_beams(params, initial_vars)
        elseif params.beam_sizer == :discrete
            params = process_discrete_beams(params)
        end
    catch e
        if e isa NoValidSectionsError
            @warn "Deflection-feasible design infeasible: $(e.msg)"
            params.minimizers = Vector{Float64}[]
        else
            rethrow(e)
        end
    end

    if params.verbose
        println("Minimums: $(params.minimums)")
    end

    # Legacy staged-bookkeeping fields — always "converged" under the
    # single-pass regime. Kept on `SlabSizingParams` for downstream API
    # compatibility.
    params.staged_converged    = true
    params.staged_n_violations = 0

    return self, params
end

"""
    optimal_beamsizer_reconcile(self, params; initial_vars=[], verbose=false)

Size beams with a lightweight outer reconciliation loop that aligns the
sizer's per-beam deflection constraint with the global-FE deflection
measured by `postprocess_slab`.

The two deflection models differ by construction:

- **Sizer constraint** (`get_deflection_constraint`, per-beam FE on a
  pin-roller span with composite `I_x` scaling) is fast and Zygote-
  differentiable.
- **Postprocess check** (`postprocess_slab`, global FE with composite
  `I_x` substituted into each beam's `Section`) captures frame coupling,
  support interaction, and load application semantics that the per-beam
  model cannot.

Empirically these diverge by a roughly constant factor across a slab
(frame-coupling tends to scale deflection uniformly). Instead of running
a global FE inside every NLP iteration, this loop:

 1. Sizes once with the current `serviceability_tighten`.
 2. Calls `postprocess_slab` and computes `max(δ / Δ_lim)`.
 3. If `max_ratio > 1 + reconcile_tol`, multiplies
    `serviceability_tighten *= max_ratio` (so the sizer effectively
    targets `L / (serviceability_lim · tighten)`), warm-starts from the
    last minimizers, and re-sizes.
 4. Stops when the ratio converges, the minimizer is empty, or after
    `params.reconcile_max_iter` iterations.

When `params.reconcile_max_iter == 0` this reduces to a single call to
`optimal_beamsizer` followed by one `postprocess_slab`, i.e. the default
single-pass regime.

# Returns

`(self, params, result, max_ratio, n_iters)` where
- `result` is the `SlabOptimResults` from the final `postprocess_slab`
  (populated for the serviceable design).
- `max_ratio` is the postprocess `δ / Δ_lim` attained at exit (≤1 means
  the sizer-postprocess gap is closed).
- `n_iters` is the number of sizer calls used (always ≥ 1).
"""
function optimal_beamsizer_reconcile(self::SlabAnalysisParams, params::SlabSizingParams;
                                     initial_vars::Vector=[], verbose::Bool=false)

    max_iter = params.reconcile_max_iter
    tol      = params.reconcile_tol
    warm     = initial_vars

    result      = nothing
    max_ratio   = 0.0
    n_iters     = 0

    for iter in 0:max_iter
        n_iters += 1
        self, params = optimal_beamsizer(self, params; initial_vars=warm)

        if isempty(params.minimizers)
            verbose && @warn "Reconciliation exited at iter $iter: no feasible sections"
            break
        end

        # Full postprocess — we need the per-beam deflections anyway to
        # decide whether to tighten. The loop exits as soon as the check
        # is satisfied.
        result = postprocess_slab(self, params)

        lim_vec = result.Δ_limit_total
        δ_vec   = result.δ_total
        ratios  = [lim_vec[i] > 0 ? δ_vec[i] / lim_vec[i] : 0.0 for i in eachindex(δ_vec)]
        max_ratio = isempty(ratios) ? 0.0 : maximum(ratios)

        verbose && println("  reconcile iter=$iter  tighten=$(round(params.serviceability_tighten, digits=3))  " *
                           "max δ/Δ=$(round(max_ratio, digits=3))  n_fail=$(count(ratios .> 1.0 + tol))")

        # Converged or disabled
        if max_ratio <= 1.0 + tol || max_iter == 0
            break
        end

        # Tighten and warm-start the next iteration from the current
        # minimizers (cheap — NLP only needs a few iterations per beam).
        params.serviceability_tighten *= max_ratio
        warm = params.minimizers
    end

    return self, params, result, max_ratio, n_iters
end

"""
    process_discrete_beams(params::SlabSizingParams) -> SlabSizingParams

Size each beam sequentially against the catalogue in `params.catalog_discrete`
using `sequential_search_sections`, with the single-limit deflection
constraint provided by `get_deflection_constraint`.
"""
function process_discrete_beams(params::SlabSizingParams)
    
    model = params.model
    beam_elements = model.elements[:beam]

    for i in 1:lastindex(beam_elements)

        println("=== $i / $(lastindex(beam_elements)) ===")
        
        n_sections = 1
        i_check = Int64[]
        beam_loads = params.load_dictionary[get_element_id(beam_elements[i])]
        beam_forces = InternalForces(beam_elements[i], beam_loads, resolution=200)

        if n_sections == 1

            # check if the section has existed before
            M_max, V_max, x_max = maximum(abs.(beam_forces.My)), maximum(abs.(beam_forces.Vy)), maximum(abs.(beam_forces.x))
            params, unique = check_uniqueness!(params, M_max, V_max, x_max)
            if !unique 
                continue 
            end

        end

        # set up variables
        fixed_variables = [true, true, true, true]
        broadcast_in, broadcast_out = generate_broadcast(fixed_variables, n_sections)

        max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))
        
        if !iszero(params.max_beam_depth) && !isinf(params.max_beam_depth)
            max_h = min(params.max_beam_depth, max_h)
        end

        # set up constraints and objective
        beam_params = get_frameoptparams(params, beam_elements[i], params.load_dictionary[get_element_id(beam_elements[i])])

        objective = generate_objective(params, beam_params, beam_forces, broadcast_out, i_check, n_sections, interpolation=akima_interpolation)
        constraints = inequality_constraints_I_symm(beam_forces, n_sections, i_check=i_check, interpolation_function=akima_interpolation, broadcast_out=broadcast_out)

        if params.deflection_limit == true

            element_loads_i = params.load_dictionary[get_element_id(beam_elements[i])]
            push!(constraints, get_deflection_constraint(beam_params, beam_forces.x[end], params,
                                                         beam_element=beam_elements[i],
                                                         element_loads=element_loads_i))

        end

        # get the sections to evaluate over
        minimum, minimum_section = sequential_search_sections(params.catalog_discrete, objective, constraints, max_h)

        # get the best section
        minimizer = get_geometry_vars(minimum_section)
        id = minimum_section.name

        println("")
        push!(params.minimizers, minimizer)
        push!(params.minimums, minimum)
        push!(params.ids, id)
        push!(params.M_maxs, M_max)
        push!(params.V_maxs, V_max)
        push!(params.x_maxs, x_max)
    end

    return params
end


"""
    process_continuous_beams(beam_elements, slab_model, minimizers, minimums, ids, M_maxs, V_maxs, x_maxs, params::SlabSizingParams, initial_vars)

Processes continuous beam sizing.
"""
function process_continuous_beams(params::SlabSizingParams, initial_vars::Vector=[])
   
    model = params.model
    beam_elements = model.elements[:beam]
   
    for i in 1:lastindex(beam_elements)
        
        println("=== $i / $(lastindex(beam_elements)) ===")

        # set up the single variable problem
        n_sections = 1
        i_check = Int64[]
        beam_loads = params.load_dictionary[get_element_id(beam_elements[i])]
        beam_forces = InternalForces(beam_elements[i], beam_loads, resolution=200)

        if n_sections == 1

            # check if the section has existed before
            M_max, V_max, x_max = maximum(abs.(beam_forces.My)), maximum(abs.(beam_forces.Vy)), maximum(abs.(beam_forces.x))
            params, unique = check_uniqueness!(params, M_max, V_max, x_max)
            if !unique
                continue
            end

        end

        # optimization
        optimization_model = Nonconvex.Model()

        # set up variables
        fixed_variables = [true, true, true, true] # h, w, tw, tf; in the case of variable optimization, these remain constant across all sections
        broadcast_in, broadcast_out = generate_broadcast(fixed_variables, n_sections)

        # -- check if there's a minimum for continuous beams
        if params.minimum_continuous == true
            min_h, min_w, min_tw, min_tf = get_geometry_vars(W_imperial("W6X8.5"))
        else
            min_h, min_w, min_tw, min_tf = [0.01, 0.01, 0.001, 0.001]
        end

        max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))
        
        if !iszero(params.max_beam_depth) && !isinf(params.max_beam_depth)
            max_h = min(params.max_beam_depth, max_h)
        end

        addvar!(optimization_model, broadcast_in([min_h, min_w, min_tw, min_tf]), broadcast_in([max_h, max_w, max_tw, max_tf])) # set max and min

        # set up constraints and objective
        beam_params = get_frameoptparams(params, beam_elements[i], params.load_dictionary[get_element_id(beam_elements[i])])
                
        objective = generate_objective(params, beam_params, beam_forces, broadcast_out, i_check, n_sections, interpolation=akima_interpolation)
        set_objective!(optimization_model, objective)
        
        constraints = inequality_constraints_I_symm(beam_forces, n_sections, i_check=i_check, interpolation_function=akima_interpolation, broadcast_out=broadcast_out)

        if params.deflection_limit == true
            element_loads_i = params.load_dictionary[get_element_id(beam_elements[i])]
            single_section_deflection = get_deflection_constraint(beam_params, beam_forces.x[end], params,
                                                                  beam_element=beam_elements[i],
                                                                  element_loads=element_loads_i)
            push!(constraints, single_section_deflection)
        end

        for constraint in constraints
            add_ineq_constraint!(optimization_model, constraint)
        end

        alg = NLoptAlg(:LD_MMA) # NLoptAlg(LD_MMA), MMA02, NLoptAlg(:LN_COBYLA), NLoptAlg(:LD_CCSAQ)
        options = NLoptOptions(
            xtol_rel=1e-2,
            xtol_abs=1e-8,
            ftol_abs=1e-2,
            ftol_rel=1e-2,
        )

        # `typeof(initial_vars) == Vector{Vector}` is Julia-type-invariant
        # (it's false for `Vector{Vector{Float64}}`), so we match the element
        # type instead to properly flow `Vector{Vector{Float64}}` warm starts.
        if isempty(initial_vars)
            init_vars = broadcast_in(get_geometry_vars(W_imperial("W6X8.5")))
        elseif eltype(initial_vars) <: AbstractString
            init_vars = broadcast_in(get_geometry_vars(W_imperial(initial_vars[i])))
        elseif eltype(initial_vars) <: AbstractVector
            init_vars = broadcast_in(initial_vars[i])
        elseif eltype(initial_vars) <: Asap.AbstractSection
            init_vars = broadcast_in(get_geometry_vars(initial_vars[i]))
        else
            init_vars = broadcast_in(get_geometry_vars(W_imperial("W6X8.5")))
        end

        # Qualify to `Nonconvex.optimize` — several loaded packages (Optim,
        # NLopt, NonconvexCore/MMA/NLopt/Ipopt/Juniper) export `optimize`.
        result = Nonconvex.optimize(optimization_model, alg, init_vars, options=options)

        println("")
        push!(params.minimizers, result.minimizer)
        push!(params.minimums, result.minimum)
        push!(params.ids, string(round(I_symm(result.minimizer...).A, digits=2)))
        push!(params.M_maxs, M_max)
        push!(params.V_maxs, V_max)
        push!(params.x_maxs, x_max)

    end

    return params

end


"""
    iterate_discrete_continuous(self::SlabAnalysisParams, geometry_dict::Dict; sections::Vector=[], max_depth::Real, deflection_limit::Bool=true, verbose::Bool=false, minimum::Bool=false, save::Bool=false)

Iterates over discrete and continuous beam sizing methods to optimize the slab design.

# Arguments
- `self::SlabAnalysisParams`: The slab analysis parameters.
- `geometry_dict::Dict`: Dictionary containing geometry data.
- `sections::Vector`: Initial sections for the beams.
- `max_depth::Real`: Maximum depth for beam sizing.
- `deflection_limit::Bool`: Whether to apply deflection limits.
- `verbose::Bool`: Whether to print detailed output.
- `minimum::Bool`: Whether to use minimum geometry constraints.
- `save::Bool`: Whether to save the results.

# Returns
- If `save` is true, returns the processed results for both discrete and continuous methods.
"""
function iterate_discrete_continuous(analysis_params::SlabAnalysisParams, sizing_params::SlabSizingParams)

    # Initialize results
    slab_results_discrete_noncollinear, slab_results_discrete_collinear, slab_results_continuous_noncollinear, slab_results_continuous_collinear = initialize_slab_results(analysis_params, sizing_params)

    initial_sections = Asap.AbstractSection[]
    initial_vars = Vector[]

    # Analyze slab
    try
        analysis_params = analyze_slab(analysis_params)            
    catch e
        if e isa AssertionError
            println("Error in slab analysis: ", e)
        else
            rethrow(e)
        end
    end

    # Optimize beam sizes discrete
    try
        sizing_params.beam_sizer = :discrete
        sizing_params = reset_SlabSizingParams(sizing_params)
        if isempty(analysis_params.slab_depths)
            return slab_results_discrete_noncollinear, slab_results_discrete_collinear, slab_results_continuous_noncollinear, slab_results_continuous_collinear
        end
        analysis_params, sizing_params = optimal_beamsizer(analysis_params, sizing_params)

        initial_sections = [to_ASAP_section(I_symm(minimizer...)) for minimizer in sizing_params.minimizers]
        initial_vars = sizing_params.minimizers
        
        slab_results_discrete_noncollinear = postprocess_slab(analysis_params, sizing_params, check_collinear=false)
        print("Noncollinear:\n")
        print_forces(slab_results_discrete_noncollinear)
        
        slab_results_discrete_collinear = postprocess_slab(analysis_params, sizing_params, check_collinear=true)
        print("Collinear:\n")
        print_forces(slab_results_discrete_collinear)
    catch e
        if isa(e, AssertionError) || isa(e, NoValidSectionsError)
            println("Error in discrete optimization: ", e)
            # Return early since we don't want to continue to continuous optimization
            return slab_results_discrete_noncollinear, slab_results_discrete_collinear, slab_results_continuous_noncollinear, slab_results_continuous_collinear
        else
            rethrow(e)
        end
    end

    # Optimize beam sizes continuous
    try
        sizing_params.beam_sizer = :continuous
        sizing_params = reset_SlabSizingParams(sizing_params)
        analysis_params, sizing_params = optimal_beamsizer(analysis_params, sizing_params, initial_vars=initial_vars)
        
        slab_results_continuous_noncollinear = postprocess_slab(analysis_params, sizing_params, check_collinear=false)
        print("Noncollinear:\n")
        print_forces(slab_results_continuous_noncollinear)
        
        slab_results_continuous_collinear = postprocess_slab(analysis_params, sizing_params, check_collinear=true)
        print("Collinear:\n")
        print_forces(slab_results_continuous_collinear)
    catch e
        if e isa AssertionError
            println("Error in discrete optimization: ", e)
        elseif isa(e, NoValidSectionsError)
            println("Caught NoValidSectionsError", e)
        else
            rethrow(e)
        end
    end

    return slab_results_discrete_noncollinear, slab_results_discrete_collinear, slab_results_continuous_noncollinear, slab_results_continuous_collinear

end
