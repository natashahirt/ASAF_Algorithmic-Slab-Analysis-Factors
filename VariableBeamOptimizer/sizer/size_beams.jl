"""
    optimal_beamsizer!(self::SlabAnalysisParams, initial_vars::Vector; max_depth::Real=21, sizing_unit::Symbol=:in, deflection_limit::Bool=true, verbose::Bool=true, minimum::Bool=true, max_assembly_depth::Bool=true)

Optimizes beam sizes based on initial variables and parameters. Returns the optimal beam sizes.
"""
function optimal_beamsizer!(self::SlabAnalysisParams, params::SlabSizingParams; initial_vars::Vector=[])
    @assert !isempty(initial_vars) "You need to input initial variables as a list of strings ['W6X8.5', ...] or vectors of floats [[4.17, 4.055, 0.27, 0.34], ...]"
    return optimal_beamsizer(self, params, initial_vars=initial_vars)
end

"""
    optimal_beamsizer(self, params; initial_vars, max_staged_iters)

Takes metric forces to size imperial beams and returns imperial beams.

When composite action is enabled and staged load columns exist in `load_df`,
an outer convergence loop verifies the sized sections against exact staged
deflection limits (L/360 live, L/240 total) using global FE solves.  Beams
that violate a limit have their minimum Ix floors tightened, and the sizer
re-runs until all beams comply or `max_staged_iters` is reached.
"""
function optimal_beamsizer(self::SlabAnalysisParams, params::SlabSizingParams;
                           initial_vars::Vector=[], max_staged_iters::Int=15)

    params.model = self.model

    if !isempty(params.M_maxs)
        self = reset_SlabAnalysisParams(self, self.model)
        params = reset_SlabSizingParams(params)
    end

    # convert model lengths to inches
    conversion_factor = convert_to_m[self.slab_units] * 1/convert_to_m[params.beam_units]
    params.area = self.area * conversion_factor^2

    # Adjust max depth based on slab depth
    if params.max_assembly_depth
        slab_depth = maximum(self.slab_depths) * conversion_factor
        params.max_beam_depth = params.max_depth - slab_depth
    end

    # Store slab depth, perimeter set, and max bay span in beam units for composite action
    if params.composite_action
        params.slab_depth_in = maximum(self.slab_depths) * conversion_factor
        params.i_perimeter = Set(self.i_perimeter)
    end
    params.max_bay_span = maximum(self.max_spans) * conversion_factor

    params.model = get_scaled_model(self, params, conversion_factor)
    params.load_dictionary = get_load_dictionary_by_id(params.model)

    # Compute collinear groups (used by MIP to enforce same-section constraint)
    if params.collinear == true
        params.collinear_groups = get_collinear_groups(params.model.elements[:beam])
    end

    has_staged = params.composite_action &&
                 params.deflection_limit &&
                 :unfactored_w_live in propertynames(params.load_df)

    # ── MIP warm-start for NLP: run discrete sizing first to seed NLP ──
    if params.beam_sizer == :continuous && isempty(initial_vars)
        try
            println("Running MIP warm-start for NLP...")
            mip_params = deepcopy(params)
            mip_params.beam_sizer = :discrete
            mip_params = process_discrete_beams_integer(mip_params)
            if !isempty(mip_params.ids)
                initial_vars = mip_params.ids
                params.mip_result = mip_params
                println("  MIP warm-start: $(length(initial_vars)) sections seeded.")
            end
        catch e
            @warn "MIP warm-start failed; NLP will use default init" exception=e
        end
    end

    # Preserve the best valid sizing result across outer-loop iterations
    best_minimizers = Vector{Float64}[]
    best_minimums   = Float64[]
    best_ids        = String[]
    best_M_maxs     = Float64[]
    best_V_maxs     = Float64[]
    best_x_maxs     = Float64[]

    # Cycle detection state
    viol_history = Int[]
    cycling      = false
    final_n_viol = 0

    for staged_iter in 1:max_staged_iters

        # ── size beams ────────────────────────────────────────────────────
        if staged_iter > 1
            params.minimizers = Vector{Float64}[]
            params.minimums   = Float64[]
            params.ids        = String[]
            params.M_maxs     = Float64[]
            params.V_maxs     = Float64[]
            params.x_maxs     = Float64[]
        end

        try
            if params.beam_sizer == :continuous
                params = process_continuous_beams_parallel(params, initial_vars=initial_vars)
            elseif params.beam_sizer == :discrete
                params = process_discrete_beams_integer(params)
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

        # Skip staged verification when not applicable
        if !has_staged
            break
        end

        # If sizing failed (e.g. infeasible MIP), restore best prior result
        if isempty(params.minimizers)
            @warn "Sizing produced no results on iteration $staged_iter; restoring best prior result."
            params.minimizers = best_minimizers
            params.minimums   = best_minimums
            params.ids        = best_ids
            params.M_maxs     = best_M_maxs
            params.V_maxs     = best_V_maxs
            params.x_maxs     = best_x_maxs
            break
        end

        # Snapshot current result as the best so far
        best_minimizers = copy(params.minimizers)
        best_minimums   = copy(params.minimums)
        best_ids        = copy(params.ids)
        best_M_maxs     = copy(params.M_maxs)
        best_V_maxs     = copy(params.V_maxs)
        best_x_maxs     = copy(params.x_maxs)

        # ── staged deflection verification (exact global FE) ──────────────
        n_viol, new_Ix_comp, new_Ix_bare = _verify_staged_deflection(params)
        final_n_viol = n_viol

        if n_viol == 0
            println("Staged deflection feasible after $staged_iter sizing iteration(s).")
            break
        end

        # ── cycle detection ───────────────────────────────────────────────
        push!(viol_history, n_viol)
        if !cycling && length(viol_history) >= 4
            recent = viol_history[end-2:end]
            if minimum(recent) >= viol_history[end-3]
                cycling = true
                println("Staged deflection cycle detected at iter $staged_iter — switching to ratchet mode.")
            end
        end

        println("Staged deflection: $n_viol beam(s) violate limits — updating Ix floors (iter $staged_iter)$(cycling ? " [ratchet]" : "")...")

        α = 0.7
        overshoot = cycling ? 1.1 : 1.0

        for (i, val) in new_Ix_comp
            target = val * overshoot
            old = get(params.min_Ix_comp, i, 0.0)
            blended = old + α * (target - old)
            params.min_Ix_comp[i] = cycling ? max(old, blended) : blended
        end
        for (i, val) in new_Ix_bare
            target = val * overshoot
            old = get(params.min_Ix_bare, i, 0.0)
            blended = old + α * (target - old)
            params.min_Ix_bare[i] = cycling ? max(old, blended) : blended
        end

        if !cycling
            # Decay floors for beams no longer in violation (disabled when cycling)
            for i in collect(keys(params.min_Ix_comp))
                if !haskey(new_Ix_comp, i)
                    params.min_Ix_comp[i] *= (1 - α)
                end
            end
            for i in collect(keys(params.min_Ix_bare))
                if !haskey(new_Ix_bare, i)
                    params.min_Ix_bare[i] *= (1 - α)
                end
            end
        end

        if staged_iter == max_staged_iters
            @warn "Staged deflection did not fully converge after $max_staged_iters iterations; $n_viol violation(s) remain."
        end
    end

    # Factored beam self-weight on the global FE model for strength demands and column axial
    if !isempty(params.minimizers)
        finalize_beam_selfweight_factored_demands!(params)
    end

    # Record convergence outcome on sizing params for downstream filtering
    params.staged_converged    = !has_staged || final_n_viol == 0
    params.staged_n_violations = final_n_viol

    return self, params
end

"""
    _verify_staged_deflection(params) -> (n_violations, Ix_comp_needed, Ix_bare_needed)

Run exact staged FE solves with the current minimizers and check L/360 / L/240.
Returns the number of violating beams and Dicts of tightened Ix requirements
for any beam that fails (keyed by beam index).
"""
function _verify_staged_deflection(params::SlabSizingParams)
    beam_elements = params.model.elements[:beam]
    n_beams = length(beam_elements)
    minimizers = params.minimizers

    E_steel = steel_ksi.E
    ρ_steel_kip = steel_ksi.ρ

    beam_ids = [get_element_id(be) for be in beam_elements]

    # Compute fresh section properties from minimizers
    sec_A   = [A_I_symm(minimizers[i]...)  for i in 1:n_beams]
    bare_Ix = [Ix_I_symm(minimizers[i]...) for i in 1:n_beams]
    sec_Iy  = [Iy_I_symm(minimizers[i]...) for i in 1:n_beams]
    sec_J   = [J_I_symm(minimizers[i]...)  for i in 1:n_beams]

    comp_Ix  = copy(bare_Ix)
    if params.composite_action && params.slab_depth_in > 0
        for i in 1:n_beams
            L_beam = beam_elements[i].length
            element_loads_i = params.load_dictionary[beam_ids[i]]
            positions_i = Float64[]
            widths_i = Float64[]
            for ld in element_loads_i
                if hasproperty(ld, :loadID)
                    row = findfirst(==(getproperty(ld, :loadID)), params.load_df.loadID)
                    if !isnothing(row)
                        push!(positions_i, ld.position)
                        push!(widths_i, params.load_df[row, :trib_width])
                    end
                end
            end
            if !isempty(widths_i)
                is_perim = i in params.i_perimeter
                comp_Ix[i] = get_I_composite_effective(
                    minimizers[i][1], minimizers[i][2], minimizers[i][3], minimizers[i][4],
                    params.slab_depth_in, E_steel, params.E_c,
                    L_beam, positions_i, widths_i; is_perimeter=is_perim)
            end
        end
    end

    # Set bare-steel sections for Solve A (fresh A, Iy, J from minimizers)
    for i in 1:n_beams
        beam_elements[i].section = Section(sec_A[i], E_steel, steel_ksi.G,
            bare_Ix[i], sec_Iy[i], sec_J[i])
    end
    update_load_values_staged!(params.model, params, load_case=:slab_dead)
    params.load_dictionary = get_load_dictionary_by_id(params.model)
    Asap.solve!(params.model, reprocess=true)

    δ_slab_dead = zeros(n_beams)
    for i in 1:n_beams
        disp = ElementDisplacements(beam_elements[i],
            params.load_dictionary[beam_ids[i]], resolution=200)
        δ_slab_dead[i] = maximum(abs.(disp.ulocal[2, :]))
    end

    # Analytical beam self-weight on bare steel
    δ_beam_dead = zeros(n_beams)
    for i in 1:n_beams
        w_sw = sec_A[i] * ρ_steel_kip
        L_in = beam_elements[i].length
        δ_beam_dead[i] = 5 * w_sw * L_in^4 / (384 * E_steel * bare_Ix[i])
    end

    # Set composite sections for Solve B (fresh A, Iy, J from minimizers)
    for i in 1:n_beams
        beam_elements[i].section = Section(sec_A[i], E_steel, steel_ksi.G,
            comp_Ix[i], sec_Iy[i], sec_J[i])
    end
    update_load_values_staged!(params.model, params, load_case=:sdl_live)
    params.load_dictionary = get_load_dictionary_by_id(params.model)
    Asap.solve!(params.model, reprocess=true)

    δ_sdl_live = zeros(n_beams)
    for i in 1:n_beams
        disp = ElementDisplacements(beam_elements[i],
            params.load_dictionary[beam_ids[i]], resolution=200)
        δ_sdl_live[i] = maximum(abs.(disp.ulocal[2, :]))
    end

    # Split SDL+LL by load ratio
    loadid_index = _build_loadid_index(params)
    δ_live = zeros(n_beams)
    δ_total = zeros(n_beams)
    for i in 1:n_beams
        w_sdl_sum = 0.0; w_live_sum = 0.0
        for ld in params.load_dictionary[beam_ids[i]]
            if hasproperty(ld, :loadID)
                row = get(loadid_index, ld.loadID, nothing)
                if !isnothing(row)
                    w_sdl_sum  += params.load_df[row, :unfactored_w_sdl]
                    w_live_sum += params.load_df[row, :unfactored_w_live]
                end
            end
        end
        w_tot = w_sdl_sum + w_live_sum
        f_live = w_tot > 0 ? w_live_sum / w_tot : 0.5
        δ_live[i] = δ_sdl_live[i] * f_live
        δ_total[i] = δ_slab_dead[i] + δ_beam_dead[i] + δ_sdl_live[i]
    end

    # ── Global deflection sanity check ────────────────────────────────
    if params.max_bay_span > 0
        max_total_δ = maximum(δ_total)
        global_lim = params.max_bay_span / 180.0
        if max_total_δ > global_lim
            @warn "Global deflection sanity: max δ_total exceeds max-bay-span/180" max_total_δ global_lim params.max_bay_span
        end
    end

    # Restore factored loads (including beam self-weight) for subsequent sizing
    refresh_factored_loads_with_beam_sw!(params)

    # Check L/360 (live) and L/240 (total) limits.
    # For violating beams, compute the composite Ix that would satisfy each
    # limit and the bare-steel Ix needed to keep dead-load deflection in check.
    # No preemptive tightening or budget splitting — the outer loop converges
    # naturally by re-solving with updated Ix floors.
    Ix_comp_needed = Dict{Int, Float64}()
    Ix_bare_needed = Dict{Int, Float64}()
    n_viol = 0

    for i in 1:n_beams
        L = beam_elements[i].length
        lim_live  = L / 360.0
        lim_total = L / 240.0
        violated = false

        # ── L/360 live ────────────────────────────────────────────────
        if lim_live > 0 && δ_live[i] > lim_live
            Ix_comp_needed[i] = comp_Ix[i] * δ_live[i] / lim_live
            violated = true
        end

        # ── L/240 total (coupled) ────────────────────────────────────
        if lim_total > 0 && δ_total[i] > lim_total
            δ_dead_portion = δ_slab_dead[i] + δ_beam_dead[i]
            residual_240 = max(lim_total - δ_dead_portion, lim_total * 0.1)

            if δ_sdl_live[i] > residual_240
                new_comp = comp_Ix[i] * δ_sdl_live[i] / residual_240
                Ix_comp_needed[i] = max(get(Ix_comp_needed, i, 0.0), new_comp)
            end

            if δ_dead_portion > lim_total
                Ix_bare_needed[i] = bare_Ix[i] * δ_dead_portion / lim_total
            end

            violated = true
        end

        if violated
            n_viol += 1
        end
    end

    return n_viol, Ix_comp_needed, Ix_bare_needed
end


"""
    process_discrete_beams!(beam_elements, slab_model, minimizers, minimums, ids, params::SlabSizingParams)

Processes discrete beam sizing.
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

        if isempty(initial_vars)
            init_vars = broadcast_in(get_geometry_vars(W_imperial("W6X8.5")))
        elseif typeof(initial_vars) == Vector{String}
            init_vars = broadcast_in(get_geometry_vars(W_imperial(initial_vars[i])))
        elseif typeof(initial_vars) == Vector{Vector}
            init_vars = broadcast_in(initial_vars[i])
        elseif typeof(initial_vars) == Vector{Section}
            init_vars = broadcast_in(get_geometry_vars(initial_vars[i]))
        end

        result = optimize(optimization_model, alg, init_vars, options = options)

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
    iterate_discrete_continuous(analysis_params, sizing_params)

Run MIP (discrete) then NLP (continuous) sizing for one slab.  The MIP
result seeds the NLP via `initial_vars`, so the MIP only runs once.

When `sizing_params.collinear == true`, same-section constraints are enforced
**during** both the MIP and NLP solves; postprocessing follows the solver flag
without overriding it.

# Returns
A 4-tuple `(discrete_noncollinear, discrete_collinear,
            continuous_noncollinear, continuous_collinear)` of `SlabOptimResults`
for backward compatibility with CSV consumers that key on `(beam_sizer, collinear)`.
When `sizing_params.collinear == true`, the noncollinear and collinear results
within each sizer are identical (solver already enforced matching sections).
"""
function iterate_discrete_continuous(analysis_params::SlabAnalysisParams, sizing_params::SlabSizingParams)

    slab_results_discrete_noncollinear, slab_results_discrete_collinear,
        slab_results_continuous_noncollinear, slab_results_continuous_collinear =
        initialize_slab_results(analysis_params, sizing_params)

    try
        analysis_params = analyze_slab(analysis_params)
    catch e
        if e isa AssertionError
            @warn "Slab analysis failed with AssertionError; returning blank results" exception=e
            return slab_results_discrete_noncollinear, slab_results_discrete_collinear,
                   slab_results_continuous_noncollinear, slab_results_continuous_collinear
        else
            rethrow(e)
        end
    end

    # ── Discrete (MIP) sizing ─────────────────────────────────────────────
    discrete_ok = false
    mip_minimizers = Vector{Float64}[]
    try
        sizing_params.beam_sizer = :discrete
        sizing_params = reset_SlabSizingParams(sizing_params)
        if isempty(analysis_params.slab_depths)
            return slab_results_discrete_noncollinear, slab_results_discrete_collinear,
                   slab_results_continuous_noncollinear, slab_results_continuous_collinear
        end
        analysis_params, sizing_params = optimal_beamsizer(analysis_params, sizing_params)

        if !isempty(sizing_params.minimizers)
            mip_minimizers = sizing_params.minimizers

            slab_results_discrete_collinear = postprocess_slab(analysis_params, sizing_params)
            slab_results_discrete_collinear.collinear = true
            print("Discrete (collinear):\n")
            print_forces(slab_results_discrete_collinear)

            slab_results_discrete_noncollinear = deepcopy(slab_results_discrete_collinear)
            slab_results_discrete_noncollinear.collinear = false

            discrete_ok = true
        end
    catch e
        if isa(e, AssertionError) || isa(e, NoValidSectionsError)
            println("Error in discrete optimization: ", e)
            return slab_results_discrete_noncollinear, slab_results_discrete_collinear,
                   slab_results_continuous_noncollinear, slab_results_continuous_collinear
        else
            rethrow(e)
        end
    end

    # ── Continuous (NLP) sizing, seeded by the MIP result above ───────────
    # Passing initial_vars skips the redundant internal MIP warm-start
    # inside optimal_beamsizer (line 62: fires only when initial_vars is empty).
    if discrete_ok
        try
            sizing_params.mip_result = deepcopy(sizing_params)
            sizing_params.beam_sizer = :continuous
            sizing_params = reset_SlabSizingParams(sizing_params)
            analysis_params, sizing_params = optimal_beamsizer(
                analysis_params, sizing_params, initial_vars=mip_minimizers)

            if !isempty(sizing_params.minimizers)
                slab_results_continuous_collinear = postprocess_slab(analysis_params, sizing_params)
                slab_results_continuous_collinear.collinear = true
                print("Continuous (collinear):\n")
                print_forces(slab_results_continuous_collinear)

                slab_results_continuous_noncollinear = deepcopy(slab_results_continuous_collinear)
                slab_results_continuous_noncollinear.collinear = false
            end
        catch e
            if e isa AssertionError
                println("Error in continuous optimization: ", e)
            elseif isa(e, NoValidSectionsError)
                println("Caught NoValidSectionsError: ", e)
            else
                rethrow(e)
            end
        end
    end

    return slab_results_discrete_noncollinear, slab_results_discrete_collinear,
           slab_results_continuous_noncollinear, slab_results_continuous_collinear
end
