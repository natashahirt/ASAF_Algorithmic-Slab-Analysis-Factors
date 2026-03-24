Zygote.@nograd I_symm

"""
    topopt_default_schedule()

Default continuation schedule for `optimize_indeterminate`.
Each stage is a named tuple with keys:
`:name, :sharpness, :penalty, :xtol_abs, :ftol_abs, :maxeval, :maxtime`.
"""
function topopt_default_schedule()
    return [
        (name="stage_1_smooth", sharpness=8.0, penalty=0.1, xtol_abs=1e-2, ftol_abs=1e-3, maxeval=1200, maxtime=180.0),
        (name="stage_2_transition", sharpness=16.0, penalty=1.0, xtol_abs=1e-3, ftol_abs=1e-4, maxeval=1600, maxtime=240.0),
        (name="stage_3_crisp", sharpness=28.0, penalty=5.0, xtol_abs=1e-4, ftol_abs=1e-5, maxeval=2200, maxtime=300.0),
    ]
end

"""
    get_indeterminate_constraint_MV(self::SlabAnalysisParams, params::SlabSizingParams; initial_vars::Vector=[], resolution::Int=200)

    Differentiably redistributes loads across beams to compute the max moment, shear, and x_max DEMAND for each beam.
    Returns a function that gives the max moment, shear, and x_max demand for each beam given a set of beam dimensions.
"""
function get_constraint_and_objective_indeterminate(self::SlabAnalysisParams, params::SlabSizingParams)

    beam_elements = params.model.elements[:beam]
    n_beams = lastindex(beam_elements)
    l = [beam.length for beam in beam_elements]
    
    conversion_factor = 1/convert_to_m[:in]

    slab_area = sum(self.load_areas)
    slab_depth = self.slab_depths[1] * conversion_factor # convert m to in
    ρ_conc_kipin3 = params.concrete_material.ρ_concrete_kipin3
    unfactored_w_applied = params.live_load + params.superimposed_dead_load
    if params.slab_dead_load > 0.0
        unfactored_w_slab = params.slab_dead_load
    elseif self.slab_type == :orth_biaxial
        unfactored_w_slab = 0.98 * ρ_conc_kipin3 + 0.02 * ρ_STEEL_KIPIN3
    else
        unfactored_w_slab = 0.99 * ρ_conc_kipin3 + 0.01 * ρ_STEEL_KIPIN3
    end
    factored_w_applied = params.live_factor * params.live_load + params.dead_factor * params.superimposed_dead_load # ksi
    factored_w_slab = params.dead_factor * unfactored_w_slab # ksi if slab_dead_load is per-area, kip/in³ otherwise (consistent with downstream use)
    unfactored_w_facade = params.façade_load # kip/in
    factored_w_facade = params.dead_factor * unfactored_w_facade # kip/in
    perimeter_indices = collect(self.i_perimeter)

    expected_load = factored_w_slab * slab_depth + factored_w_applied
    expected_total_load = expected_load * slab_area

    println("Expected load: $(round(expected_load, digits=6)) ksi")
    println("Expected total load: $(round(expected_total_load, digits=2)) kips")

    cell_data = Dict{Int, CellData}()

    # make the cell data dictionary
    for row in eachrow(self.raster_df)
        centerpoint_id = row.centerpoint_id
        if haskey(cell_data, centerpoint_id)
            # Accumulate data for same centerpoint
            existing = cell_data[centerpoint_id]
            cell_data[centerpoint_id] = CellData(
                existing.area + row.areas,
                vcat(existing.distances, row.distance),
                vcat(existing.beam_idxs, row.beam_idx),
                vcat(existing.beam_t, row.t)
            )
        else
            cell_data[centerpoint_id] = CellData(row.areas, [row.distance], [row.beam_idx], [row.t])
        end
    end

    # Continuation knobs (updated stage-by-stage by optimize_indeterminate).
    suppress_sharpness = Ref(20.0)
    penalty_weight = Ref(1.0)
    δ_allow = l ./ params.serviceability_lim ./ params.deflection_reduction_factor

    function evaluate_design_response(x::AbstractVector)
        x_vec = collect(x)

        h = suppress_if_small.(x_vec[1:n_beams], :h; sharpness=suppress_sharpness[])
        w = suppress_if_small.(x_vec[n_beams+1:2n_beams], :w; sharpness=suppress_sharpness[])
        tw = suppress_if_small.(x_vec[2n_beams+1:3n_beams], :tw; sharpness=suppress_sharpness[])
        tf = suppress_if_small.(x_vec[3n_beams+1:4n_beams], :tf; sharpness=suppress_sharpness[])

        max_loads_strength = solve_indeterminate_system(h, w, tw, tf, params.model, cell_data, factored_w_slab, factored_w_applied, slab_depth,
                                                        perimeter_indices=perimeter_indices, facade_line_load=factored_w_facade)
        max_Mu = max_loads_strength[1:n_beams]
        max_Vu = max_loads_strength[n_beams+1:2*n_beams]

        max_loads_service = solve_indeterminate_system(h, w, tw, tf, params.model, cell_data, unfactored_w_slab, unfactored_w_applied, slab_depth,
                                                       perimeter_indices=perimeter_indices, facade_line_load=unfactored_w_facade)
        max_δ = max_loads_service[2n_beams+1:3n_beams]

        section_list = [I_symm(h[i], w[i], tw[i], tf[i]) for i in 1:n_beams]
        A = getfield.(section_list, :A)
        Mn = getfield.(section_list, :Mn)
        Vn = getfield.(section_list, :Vn)

        return (h=h, w=w, tw=tw, tf=tf, max_Mu=max_Mu, max_Vu=max_Vu, max_δ=max_δ, A=A, Mn=Mn, Vn=Vn)
    end

    # Return the max moment, shear, and x_max DEMAND for each beam
    function indeterminate_constraint(x)
        eval = evaluate_design_response(x)
        max_Mu, max_Vu, max_δ = eval.max_Mu, eval.max_Vu, eval.max_δ
        Mn, Vn = eval.Mn, eval.Vn

        if params.verbose
            println("h: $(eval.h[1]), w: $(eval.w[1]), tw: $(eval.tw[1]), tf: $(eval.tf[1]), A: $(eval.A[1]), Mn: $(eval.Mn[1]), Vn: $(eval.Vn[1])")
        end

        ϕ_b = 0.9 # resistance factor for moment, Mu ≤ ϕ_b * Mn
        ϕ_v = 0.9 # resistance factor for shear, Vu ≤ ϕ_v * Vn

        """ start with valid geometry """

        """# ensure area is positive
        A = .-A

        # make sure 2 * tf ≤ h 
        Δ_height = (2 * tf .- h)   # Normalized height difference relative to total height
        # ratio of web area to compression flange area ≤ 10
        h_w = h .- 2 * tf
        A_w = tw .* h_w
        A_fc = tf .* w

        # Geometric constraints
        Δ_height_constraint = Δ_height  # If Δ_height should be ≤ 0
        ratio_w_to_f_constraint = A_w ./ A_fc .- 10  # Should be ≤ 0
        ratio_h_to_tw_constraint = h ./ tw .- 260    # Should be ≤ 0"""

        # Force constraints
        moment_constraint = max_Mu .- (ϕ_b * Mn)     # Should be ≤ 0
        shear_constraint = max_Vu .- (ϕ_v * Vn)      # Should be ≤ 0

        # Check deflection limit per beam (L / serviceability_lim), using local vertical deflection.
        deflection_constraint = max_δ .- δ_allow  # Should be ≤ 0
        if params.verbose
            println("Max beam deflection: $(maximum(max_δ)) inches")
        end

        # Return as a vector
        return vcat(
            moment_constraint,
            shear_constraint,
            deflection_constraint
        )

    end

    function indeterminate_objective(x)
        eval = evaluate_design_response(x)
        max_Mu, max_Vu, max_δ = eval.max_Mu, eval.max_Vu, eval.max_δ
        A, Mn, Vn = eval.A, eval.Mn, eval.Vn

        ϕ_b = 0.9
        ϕ_v = 0.9

        # Penalty for moment violations
        moment_violations = max.(0, max_Mu .- (ϕ_b .* Mn))
        moment_penalty = sum(moment_violations.^2)  # or ^p for p>1

        # Penalty for shear violations
        shear_violations = max.(0, max_Vu .- (ϕ_v .* Vn))
        shear_penalty = sum(shear_violations.^2)

        # Penalty for deflection violation (per beam span limit)
        deflection_violation = max.(0, max_δ .- δ_allow)
        deflection_penalty = sum(deflection_violation.^2)

        volume = get_volume(A, beam_elements)

        λ = penalty_weight[]
        objective = volume + λ * (moment_penalty + shear_penalty + deflection_penalty)

        if params.verbose
            println("Volume: $volume, Penalty: $(moment_penalty + shear_penalty + deflection_penalty)")
        end

        return objective

    end 

    function set_continuation_params!(; sharpness::Union{Nothing,Float64}=nothing, penalty::Union{Nothing,Float64}=nothing)
        if !isnothing(sharpness)
            suppress_sharpness[] = sharpness
        end
        if !isnothing(penalty)
            penalty_weight[] = penalty
        end
        return nothing
    end

    function get_load_dictionary(params, minimizers)

        params.model.loads = Asap.AbstractLoad[]

        h = minimizers[1:n_beams]
        w = minimizers[n_beams+1:2n_beams] 
        tw = minimizers[2n_beams+1:3n_beams]
        tf = minimizers[3n_beams+1:4n_beams]

        beam_ids = [get_element_id(params.model.elements[:beam][i]) for i in 1:n_beams]
        stiffness_beams = compute_beam_stiffnesses(h, w, tw, tf, l)
        triplets = collect_load_triplets(cell_data, stiffness_beams, factored_w_slab, factored_w_applied, slab_depth)
        load_dict = Dict{Tuple{Int, Int}, Vector{Asap.AbstractLoad}}()
        params.load_df = DataFrame(
            loadID = Int[],
            area = Float64[],
            volume = Float64[],
            depth = Float64[],
            factored_w_slab = Float64[],
            unfactored_w_slab = Float64[],
            factored_w_applied = Float64[],
            unfactored_w_applied = Float64[],
        )

        next_load_id = 1
        q_factored_total = factored_w_slab * slab_depth + factored_w_applied
        q_unfactored_total = unfactored_w_slab * slab_depth + unfactored_w_applied

        for (i, t, w) in triplets
            load = Asap.PointLoad(beam_elements[i], t, [0,0,-w])
            load.loadID = next_load_id
            push!(params.model.loads, load)
            if !haskey(load_dict, beam_ids[i])
                load_dict[beam_ids[i]] = Asap.AbstractLoad[]
            end
            push!(load_dict[beam_ids[i]], load)

            factored_slab_component = q_factored_total > 0 ? w * (factored_w_slab * slab_depth) / q_factored_total : 0.0
            factored_applied_component = q_factored_total > 0 ? w * factored_w_applied / q_factored_total : 0.0
            unfactored_slab_component = q_unfactored_total > 0 ? w * (unfactored_w_slab * slab_depth) / q_unfactored_total : 0.0
            unfactored_applied_component = q_unfactored_total > 0 ? w * unfactored_w_applied / q_unfactored_total : 0.0
            area_equiv = q_factored_total > 0 ? w / q_factored_total : 0.0

            push!(params.load_df, (
                next_load_id,
                area_equiv,
                area_equiv * slab_depth,
                slab_depth,
                factored_slab_component,
                unfactored_slab_component,
                factored_applied_component,
                unfactored_applied_component,
            ))
            next_load_id += 1
        end

        # Add perimeter facade line loads to final model/load dictionary.
        if !isempty(perimeter_indices) && !iszero(factored_w_facade)
            for i in perimeter_indices
                line_load = LineLoad(beam_elements[i], [0,0,-factored_w_facade])
                push!(params.model.loads, line_load)
                if !haskey(load_dict, beam_ids[i])
                    load_dict[beam_ids[i]] = Asap.AbstractLoad[]
                end
                push!(load_dict[beam_ids[i]], line_load)
            end
        end

        # Sum up the loads from the triplets
        total_load = sum(w for (_, _, w) in triplets)
        println("Final load from triplets: $total_load")

        # Sum up the loads from the load_dict
        total_load_dict = sum(sum(load.value[3] for load in loads) for (_, loads) in load_dict)
        println("Final load from load_dict: $total_load_dict")

        params.load_dictionary = load_dict

        return params

    end

    return indeterminate_constraint, indeterminate_objective, get_load_dictionary, set_continuation_params!

end

function optimize_indeterminate(self::SlabAnalysisParams, params::SlabSizingParams;
                                initial_vars::Vector=[],
                                continuation_schedule::Union{Nothing,Vector}=nothing)

    # Initialize
    n_beams = length(self.model.elements[:beam])

    if !isempty(params.M_maxs)
        self = reset_SlabAnalysisParams(self, self.model)
        params = reset_SlabSizingParams(params)
    end

    self.load_type = :indeterminate
    self.slab_type = :isotropic
    self = analyze_slab(self)

    conversion_factor = 1/convert_to_m[:in]
    params.area = self.area * conversion_factor^2
    self.areas .*= conversion_factor^2
    self.load_areas .*= conversion_factor^2

    params.model = get_scaled_model(self, params, conversion_factor)
    params.load_dictionary = get_load_dictionary_by_id(params.model)

    # Give it a warm start
    """initial_vars = [get_geometry_vars(W_imperial("W8X35")) for _ in 1:length(self.model.elements[:beam])]
    self, beam_sizing_params = optimal_beamsizer(self, params, initial_vars = initial_vars) # differentiable, uses Ipopt
    initial_vars = beam_sizing_params.minimizers"""

    if isempty(initial_vars)
        initial_vars = [get_geometry_vars(W_imperial("W6X8.5")) for _ in 1:length(self.model.elements[:beam])]
    end
    
    # Set up optimization
    optimization_model = Nonconvex.Model()

    # Add design variables
    min_h, min_w, min_tw, min_tf = [0.01, 0.01, 0.001, 0.001]

    max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))
    max_h = minimum([params.max_depth, max_h])

    lower_bounds = vcat(fill(min_h, n_beams), fill(min_w, n_beams), fill(min_tw, n_beams), fill(min_tf, n_beams))
    upper_bounds = vcat(fill(max_h, n_beams), fill(max_w, n_beams), fill(max_tw, n_beams), fill(max_tf, n_beams))

    # Construct x0 from initial_vars
    h = [vars[1] for vars in initial_vars]
    w = [vars[2] for vars in initial_vars]
    tw = [vars[3] for vars in initial_vars]
    tf = [vars[4] for vars in initial_vars]
    x0 = vcat(h, w, tw, tf)

    h_vars = [Nonconvex.addvar!(optimization_model, min_h, max_h, init=initial_vars[i][1]) for i in 1:n_beams]
    w_vars = [Nonconvex.addvar!(optimization_model, min_w, max_w, init=initial_vars[i][2]) for i in 1:n_beams]
    tw_vars = [Nonconvex.addvar!(optimization_model, min_tw, max_tw, init=initial_vars[i][3]) for i in 1:n_beams]
    tf_vars = [Nonconvex.addvar!(optimization_model, min_tf, max_tf, init=initial_vars[i][4]) for i in 1:n_beams]    

    constraint, objective, get_load_dictionary, set_continuation_params! = get_constraint_and_objective_indeterminate(self, params)

    Nonconvex.set_objective!(optimization_model, objective)
    Nonconvex.add_ineq_constraint!(optimization_model, constraint)

    area_suppress_sharpness = Ref(20.0)

    # Keep all beams above a minimal physical section.
    min_area_threshold = 0.383 # ~1 order of magnitude smaller than the smallest beam

    for i in 1:n_beams
        Nonconvex.add_ineq_constraint!(optimization_model, x -> begin
            h = suppress_if_small(x[i], :h; sharpness=area_suppress_sharpness[])
            w = suppress_if_small(x[i + n_beams], :w; sharpness=area_suppress_sharpness[])
            tw = suppress_if_small(x[i + 2n_beams], :tw; sharpness=area_suppress_sharpness[])
            tf = suppress_if_small(x[i + 3n_beams], :tf; sharpness=area_suppress_sharpness[])
            
            # Calculate the area
            A = A_I_symm(h, w, tw, tf)
            
            # Enforce A >= min_area_threshold
            return min_area_threshold - A
        end)
    end

    # Continuation schedule: smooth -> transition -> crisp
    continuation_schedule = isnothing(continuation_schedule) ? topopt_default_schedule() : continuation_schedule
    required_stage_keys = (:name, :sharpness, :penalty, :xtol_abs, :ftol_abs, :maxeval, :maxtime)
    if isempty(continuation_schedule)
        throw(ArgumentError("continuation_schedule must contain at least one stage"))
    end
    for (k, stage) in enumerate(continuation_schedule)
        for field in required_stage_keys
            if !(field in keys(stage))
                throw(ArgumentError("continuation_schedule[$k] is missing key :$field"))
            end
        end
    end

    local_alg = NLoptAlg(:LD_MMA)
    fallback_alg = NLoptAlg(:LN_COBYLA)
    optimization_result = nothing
    x_current = copy(x0)
    for stage in continuation_schedule
        println("Running $(stage.name) with sharpness=$(stage.sharpness), penalty=$(stage.penalty)")
        set_continuation_params!(sharpness=stage.sharpness, penalty=stage.penalty)
        area_suppress_sharpness[] = stage.sharpness
        local_options = NLoptOptions(
            xtol_abs = stage.xtol_abs,
            ftol_abs = stage.ftol_abs,
            maxeval = stage.maxeval,
            maxtime = stage.maxtime
        )
        try
            optimization_result = Nonconvex.optimize(optimization_model, local_alg, x_current, options=local_options)
        catch e
            if e isa MethodError
                @warn "Stage $(stage.name) failed with gradient-based MMA; falling back to COBYLA" exception=(e, catch_backtrace())
                optimization_result = Nonconvex.optimize(optimization_model, fallback_alg, x_current, options=local_options)
            else
                rethrow(e)
            end
        end
        x_current = copy(optimization_result.minimizer)
    end

    final_sharpness = last(continuation_schedule).sharpness
    h = suppress_if_small.(optimization_result.minimizer[1:n_beams], :h; sharpness=final_sharpness)
    w = suppress_if_small.(optimization_result.minimizer[n_beams+1:2n_beams], :w; sharpness=final_sharpness)
    tw = suppress_if_small.(optimization_result.minimizer[2n_beams+1:3n_beams], :tw; sharpness=final_sharpness)
    tf = suppress_if_small.(optimization_result.minimizer[3n_beams+1:4n_beams], :tf; sharpness=final_sharpness)

    # postprocess the results
    minimizers = [[h[i], w[i], tw[i], tf[i]] for i in 1:n_beams]

    params.minimizers = minimizers
    areas = [I_symm(minimizers[i]...).A for i in 1:n_beams]
    params.minimums = [areas[i] * params.model.elements[:beam][i].length * 39.3701 for i in 1:n_beams]
    params.ids = string.(round.(areas, digits=2))

    params = get_load_dictionary(params, optimization_result.minimizer)

    return minimizers, params

end