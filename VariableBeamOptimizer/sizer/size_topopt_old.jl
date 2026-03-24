mutable struct BeamLoadGroup
    l::Float64
    x::Vector{Float64}
    t::Vector{Float64}
    areas::Vector{Float64}
    volumes::Vector{Float64}
    distances::Vector{Float64}
    Flocal::Vector{Float64}
    moment_function::Function
    shear_function::Function
    xinc::Vector{Float64}
    total_area::Float64
    total_volume::Float64
end

function process_continuous_beams_topopt_upsidedown(self::SlabAnalysisParams, params::SlabSizingParams; initial_vars::Vector=[], resolution::Int=200)

    self.load_type = :indeterminate
    self.slab_type = :isotropic
    self = analyze_slab(self)

    model = self.model
    beam_elements = model.elements[:beam]
    n_beams = length(beam_elements)

    # set up constants
    
    factored_w_applied = params.live_factor * params.live_load + params.dead_factor * params.superimposed_dead_load # ksi
    factored_w_slab = params.dead_factor * (0.99 * ρ_CONCRETE_KIPIN3 + 0.01 * ρ_STEEL_KIPIN3) # kilopounds per cubic inch, concrete + reinforcing steel

    ϕ_b = ϕ_v = 0.9
    Fy = steel_ksi.Fy

    # set up beam dataframe (may be obsolete once we have function-based moments and shears)    
    beams_df = DataFrame(i=Int64[], beamloadgroup=BeamLoadGroup[])

    # get the internal forces for each beam
    for i in 1:lastindex(beam_elements)
        beam_id = get_element_id(beam_elements[i])
        row_indices = findall(x -> x == beam_id, self.raster_df.beam_id)
        uglobal = [beam_elements[i].nodeStart.displacement; beam_elements[i].nodeEnd.displacement]
        release = AsapToolkit.get_release_type(beam_elements[i])
        r2dof = AsapToolkit.release2DOF[release]
        Flocal = (beam_elements[i].R * beam_elements[i].K * uglobal) .* r2dof
        xinc = collect(range(0, beam_elements[i].length, resolution))
        total_area = sum(self.raster_df.areas[row_indices])
        total_volume = sum(self.raster_df.volumes[row_indices])

        moment_function = AsapToolkit.MPointLoad[release]
        shear_function = AsapToolkit.VPointLoad[release]
        beamloadgroup = BeamLoadGroup(  beam_elements[i].length, 
                                        self.raster_df.x[row_indices], 
                                        self.raster_df.t[row_indices], 
                                        self.raster_df.areas[row_indices], 
                                        self.raster_df.volumes[row_indices], 
                                        self.raster_df.distance[row_indices], 
                                        Flocal,
                                        moment_function,
                                        shear_function,
                                        xinc,
                                        total_area,
                                        total_volume)
        push!(beams_df, [i, beamloadgroup])
    end

    # set up jump model

    jump_model = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(jump_model, "print_level", 1)

    # define geometry functions
    JuMP.register(jump_model, :A_I_asymm, 6, A_I_asymm, autodiff=true)
    JuMP.register(jump_model, :Iy_I_asymm, 6, Iy_I_asymm, autodiff=true)
    JuMP.register(jump_model, :parallel_axis_I, 3, parallel_axis_I, autodiff=true)
    JuMP.register(jump_model, :Zx_I_asymm, 6, Zx_I_asymm, autodiff=true)

    # define constraint functions
    function constraint_area(h, w, tw, tf)
        return -A_I_asymm(h, w, w, tw, tf, tf)  # negative area constraint
    end

    function constraint_height(h, tf)
        return 2 * tf - h  # ensure 2*tf ≤ h
    end

    function constraint_web_flange_ratio(h, w, tw, tf)
        h_w = h - 2 * tf
        A_w = tw * h_w
        A_fc = tf * w
        return A_w / A_fc - 10  # ratio ≤ 10
    end

    function constraint_web_slenderness(h, tw)
        return h / tw - 260  # h/tw ≤ 260
    end

    function constraint_flange_ratio_min(h, w, tw, tf)
        # Calculate I_yc and Iy directly
        I_yc = parallel_axis_I(w, tf, 0)
        Iy = Iy_I_asymm(h, w, w, tw, tf, tf)
        return 0.1 - I_yc/Iy  # lower bound constraint
    end

    function constraint_flange_ratio_max(h, w, tw, tf)
        I_yc = parallel_axis_I(w, tf, 0)
        Iy = Iy_I_asymm(h, w, w, tw, tf, tf)
        return I_yc/Iy - 0.9  # upper bound constraint
    end

    function Vn_I_symm(h, tw, Fy)
        Aw = h * tw
        Vn = 0.6 * Fy * Aw
        return Vn
    end

    function Ix_I_symm(h, w, tw, tf)
        h_w = h - 2*tf
        Af = w * tf
        Aw = h_w * tw
        A = 2Af + Aw

        y_b = tf / 2
        y_w = tf + h_w / 2
        y_t = h - tf / 2

        y_c = ((Af * y_b) + (Aw * y_w) + (Af * y_t)) / A

        I_b = parallel_axis_I(w, tf, (y_b - y_c))
        I_w = parallel_axis_I(tw, h_w, (y_w - y_c))
        I_t = parallel_axis_I(w, tf, (y_t - y_c))

        Ix = I_b + I_w + I_t
        return Ix
    end

    function update_raster_loads(h, w, tw, tf, i)
        group = beams_df.beamloadgroup[i]
        
        Ix = Ix_I_symm(h, w, tw, tf)
        l = group.l
        
        # Extract and transform forces
        _, Vystart, Mystart, _, _ = group.Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]
        
        My = Vector{Any}(undef, length(group.xinc))
        Vy = Vector{Any}(undef, length(group.xinc))
        My .= Vystart .* group.xinc .- Mystart
        Vy .= Vystart

        # Precompute stiffnesses and probabilities
        stiffnesses = [distance == 0 ? 1 : 1 / distance^3 + 1 / (Ix / l^3) for distance in group.distances]
        total_stiffness = sum(stiffnesses)
        probabilities = stiffnesses ./ total_stiffness        

        # Calculate areas, volumes, and loads
        areas = probabilities .* group.areas
        volumes = probabilities .* group.volumes
        loads = factored_w_slab * volumes + factored_w_applied * areas

        # Update moments and shears
        for (i, x) in enumerate(group.x)
            My .= My .+ group.moment_function.(loads[i], l, group.xinc, group.t[i])
            Vy .= Vy .+ group.shear_function.(loads[i], l, group.xinc, group.t[i])
        end

        M_max = maximum(abs.(My))
        V_max = maximum(abs.(Vy))

        println("M_max: $M_max, V_max: $V_max")
        
        return M_max, V_max
    end

    # Register all functions needed for constraints
    JuMP.register(jump_model, :constraint_area, 4, constraint_area, autodiff=true)
    JuMP.register(jump_model, :constraint_height, 2, constraint_height, autodiff=true)
    JuMP.register(jump_model, :constraint_web_flange_ratio, 4, constraint_web_flange_ratio, autodiff=true)
    JuMP.register(jump_model, :constraint_web_slenderness, 2, constraint_web_slenderness, autodiff=true)
    JuMP.register(jump_model, :constraint_flange_ratio_min, 4, constraint_flange_ratio_min, autodiff=true)
    JuMP.register(jump_model, :constraint_flange_ratio_max, 4, constraint_flange_ratio_max, autodiff=true)
    JuMP.register(jump_model, :Vn_I_symm, 3, Vn_I_symm, autodiff=true)
    JuMP.register(jump_model, :Ix_I_symm, 4, Ix_I_symm, autodiff=true)
    JuMP.register(jump_model, :update_raster_loads, 5, update_raster_loads, autodiff=true)

    # SET UP THE MODEL

    # jump: variable values
    if params.minimum_continuous == true
        min_h, min_w, min_tw, min_tf = get_geometry_vars(W_imperial("W6X8.5"))
    else
        min_h, min_w, min_tw, min_tf = [0.01, 0.01, 0.001, 0.001]
    end

    max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))

    if !iszero(params.max_beam_depth) && !isinf(params.max_beam_depth)
        max_h = min(params.max_beam_depth, max_h)
    end

    # jump: initialize variables IF available

    init_vars_array = Vector{Vector{Float64}}(undef, lastindex(beam_elements))
    for i in 1:lastindex(beam_elements)
        if isempty(initial_vars)
            init_vars_array[i] = get_geometry_vars(W_imperial("W6X8.5"))
        elseif typeof(initial_vars) == Vector{String}
            init_vars_array[i] = get_geometry_vars(W_imperial(initial_vars[i]))
        elseif typeof(initial_vars) == Vector{Vector}
            init_vars_array[i] = initial_vars[i]
        elseif typeof(initial_vars) == Vector{Section}
            init_vars_array[i] = get_geometry_vars(initial_vars[i])
        end
    end

    # Create variables with beam-specific initial values
    # each variable is a 1-dimensional vector of variables
    JuMP.@variable(jump_model, min_h <= h[i=1:lastindex(beam_elements)] <= max_h, start=init_vars_array[i][1])
    JuMP.@variable(jump_model, min_w <= w[i=1:lastindex(beam_elements)] <= max_w, start=init_vars_array[i][2])
    JuMP.@variable(jump_model, min_tw <= tw[i=1:lastindex(beam_elements)] <= max_tw, start=init_vars_array[i][3])
    JuMP.@variable(jump_model, min_tf <= tf[i=1:lastindex(beam_elements)] <= max_tf, start=init_vars_array[i][4])

    # Define objective function using the separate h, w, tw, tf variables
    JuMP.@NLobjective(jump_model, Min, sum(
        A_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) * beams_df.beamloadgroup[i].l for i in 1:n_beams)
    )

    # Add constraints using direct mathematical expressions
    for i in 1:n_beams
        # Geometric constraints
        JuMP.@NLconstraint(jump_model, constraint_area(h[i], w[i], tw[i], tf[i]) <= 0)
        JuMP.@NLconstraint(jump_model, constraint_height(h[i], tf[i]) <= 0)
        JuMP.@NLconstraint(jump_model, constraint_web_flange_ratio(h[i], w[i], tw[i], tf[i]) <= 0)
        JuMP.@NLconstraint(jump_model, constraint_web_slenderness(h[i], tw[i]) <= 0)
        
        # Flange ratio constraints (returns two constraints)
        #JuMP.@NLconstraint(jump_model, constraint_flange_ratio_min(h[i], w[i], tw[i], tf[i]) <= 0)
        #JuMP.@NLconstraint(jump_model, constraint_flange_ratio_max(h[i], w[i], tw[i], tf[i]) <= 0)

        # Strength constraints
        # Note: Using simplified versions since we can't use the full I_symm struct
        JuMP.@NLconstraint(jump_model, 
            update_raster_loads(h[i], w[i], tw[i], tf[i], i)[1] - ϕ_b * Fy * Zx_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)  # Moment capacity
        JuMP.@NLconstraint(jump_model, 
            update_raster_loads(h[i], w[i], tw[i], tf[i], i)[2] - ϕ_v * Vn_I_symm(h[i], tw[i], Fy) <= 0)  # Shear capacity
    end

    JuMP.optimize!(jump_model)

    # Extract optimized beam sizes
    optimal_h = JuMP.value.(h)
    optimal_w = JuMP.value.(w)
    optimal_tw = JuMP.value.(tw)
    optimal_tf = JuMP.value.(tf)

    minimizer = [[optimal_h[i], optimal_w[i], optimal_tw[i], optimal_tf[i]] for i in 1:n_beams]
    minimum = [A_I_symm(minimizer[i]...) * beams_df.beamloadgroup[i].l for i in 1:n_beams]

    println("Sized $(length(minimizer)) beams using simultaneous optimization.")

    params.minimizers = minimizer
    params.minimums = minimum
    params.ids = string.(round.(minimum, digits=2))
    params.M_maxs = beams_df.M_max
    params.V_maxs = beams_df.V_max
    params.x_maxs = beams_df.x_max

    return params

end