using Zygote: @nograd

function get_internal_forces_slab(slab_params::SlabAnalysisParams, beam_sizing_params::SlabSizingParams)

    slab_params = analyze_slab(slab_params, area = 100);
    conversion_factor = convert_to_m[slab_params.slab_units] * 1/convert_to_m[beam_sizing_params.beam_units]
    slab_params.area = slab_params.area * conversion_factor^2

    slab_params.model = get_scaled_model(slab_params, beam_sizing_params, conversion_factor);
    slab_params.load_dictionary = get_load_dictionary_by_id(slab_params.model);

    # Preallocate arrays for beam data
    beam_elements = slab_params.model.elements[:beam]
    n_beams = length(beam_elements)
    beam_ids = Vector{Any}(undef, n_beams)
    M_maxs = []
    V_maxs = []
    x_maxs = []

    # Get the internal forces for each beam
    for i in 1:n_beams
        beam_ids[i] = get_element_id(beam_elements[i])
        beam_loads = slab_params.load_dictionary[beam_ids[i]]
        M_max, V_max, x_max = internalforces_M_V_all(beam_elements[i], beam_loads, resolution=200)
        push!(M_maxs, M_max)
        push!(V_maxs, V_max)
        push!(x_maxs, x_max)
    end

    return M_maxs, V_maxs, x_maxs
end


"""
    get_min_height(moment::Real, shear::Real, width::Real; type::Symbol=:C16)

Calculate the minimum required height for a lumber beam based on bending and shear requirements.
"""
function get_min_height(moment::Real, shear::Real, width::Real; type::Symbol=:C16)
    material = lumber_dict[type]
    
    # Bending governs: σ = M*c / I = 6M / (b*h²), solve for h
    h_min_Mu = sqrt(6 * moment / (material.Fx * width)) # maximum fiber bending stress
    # Shear governs: τ_avg = 1.5*V / (b*h), solve for h
    h_min_Vu = 1.5 * shear / (material.Fxy * width) # max shear stress = 1.5 * average shear stress in rectangular cross sections

    return max(h_min_Mu, h_min_Vu)
end

function get_min_heights(M_maxs::Vector, V_maxs::Vector; b::Real=2.0, type::Symbol=:C24)
    min_heights = []
    for i in 1:length(M_maxs)
        push!(min_heights, get_min_height.(M_maxs[i], V_maxs[i], b, type=type))
    end
    return min_heights
end

"""
    fit_segment(x_points::Vector{<:Real}, y_points::Vector{<:Real}, x_knots::Vector{<:Real}; 
               moment::Vector{<:Real}=Float64[], b::Real=2.0, type::Symbol=:C24)

Fit a piecewise-linear function to the given points with bending stress constraints.
"""
function fit_segments(x_points::Vector{<:Real}, y_points::Vector{<:Real}, x_knots::Vector{<:Real}; moment::Vector{<:Real}=Float64[], b::Real=2.0, type::Symbol=:C24)

    x_knots = unique(sort(map(x -> min(x, x_points[end]), x_knots)))
    K, n = length(x_knots), length(x_points)
    seg_id = clamp.([findlast(x -> x <= xi, x_knots) for xi in x_points], 1, K - 1)
    segment_indices = isempty(moment) ? [] : [findall(x -> x_knots[j] ≤ x ≤ x_knots[j+1], x_points) for j in 1:(K-1)]

    function solve_model(use_nonlinear::Bool)
        model = use_nonlinear ? JuMP.Model(Ipopt.Optimizer) : JuMP.Model(Gurobi.Optimizer)
        JuMP.set_silent(model)
        JuMP.@variable(model, y_k[1:K])

        # Upper envelope constraints
        for i in 1:n
            j = seg_id[i]
            xj, xj1 = x_knots[j], x_knots[j+1]
            α = (x_points[i] - xj) / (xj1 - xj)
            JuMP.@constraint(model, (1 - α) * y_k[j] + α * y_k[j+1] ≥ y_points[i])
        end

        # Nonlinear slope constraints only if needed
        if use_nonlinear && !isempty(moment)
            Fy = lumber_dict[type].Fy
            for j in 1:K-1
                xj, xj1 = x_knots[j], x_knots[j+1]
                dx = xj1 - xj
                M_seg = maximum(moment[segment_indices[j]])
                denom = 12 * M_seg
                JuMP.@NLconstraint(model,
                ((y_k[j+1] - y_k[j]) / dx)^2 ≤
                ((Fy * 2 * b * ((y_k[j+1] + y_k[j]) / 2)^2) + 1e-8) / denom
                )
            end
        end

        JuMP.@objective(model, Min, sum(y_k))
        JuMP.optimize!(model)

        return model, JuMP.value.(y_k)
    end

    # First pass: linear only
    model, y_vals = solve_model(false)

    # Check slope violations if moment data is given
    if !isempty(moment)
        Fy = lumber_dict[type].Fy
        slope_violation = false

        for j in 1:K-1
            dx = x_knots[j+1] - x_knots[j]
            dy = y_vals[j+1] - y_vals[j]
            
            # Skip if no moment data for this segment
            if isempty(segment_indices[j])
                continue
            end
            
            M_seg = maximum(moment[segment_indices[j]])
            denom = 12 * M_seg

            s_max = sqrt(((Fy * 2 * b * ((y_vals[j] + y_vals[j+1]) / 2)^2) + 1e-8) / denom)
            slope_violation |= abs(dy / dx) > s_max
        end

        if slope_violation
            # Re-solve with nonlinear constraints
            model, y_vals = solve_model(true)
        end
    end

    # Final check for failure
    if JuMP.termination_status(model) != JuMP.MOI.OPTIMAL
        return x_knots, fill(Inf, K), fill(0.0, K - 1)
    end

    return x_knots, y_vals, diff(x_knots)
    
end
