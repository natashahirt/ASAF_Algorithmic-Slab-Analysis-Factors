# everything that goes into this function should be in inches and kips
function get_rc_objective(zones::Vector{LongitudinalZone}, material::ReinforcedConcrete; max_h::Float64=20., ϕb::Float64 = 0.9, ϕv::Float64 = 0.75, d′=1., clear_side=1., clear_bottom=2.5, clear_top=1., α=π/4, use_critical_Vn::Bool=true, rounded=false, plot::Bool=false)

    # Get beam length
    beam_length = zones[end].x_range[end]

    # Approximate shear diameter
    Av_default = REBAR_TYPES[4].diameter

    d′ = d′ + Av_default
    current_As = current_As′ = 0

    function get_rebar_layout_tension_optim(b, zones)
        max_bar_size = 11
        bars_per_layer = 6
        min_bar_number = 2
        max_unique_bars = 3
        max_layers = 3
        compression=false
        return get_rebar_layout(b, zones, max_layers, compression, d′, ϕb, clear_side, max_bar_size, bars_per_layer, min_bar_number, max_unique_bars)
    end

    function get_rebar_layout_compression_optim(b, zones)
        max_bar_size = 11
        bars_per_layer = 6
        min_bar_number = 2
        max_unique_bars = 1
        max_layers = 1
        compression=true
        return get_rebar_layout(b, zones, max_layers, compression, d′, ϕb, clear_side, max_bar_size, bars_per_layer, min_bar_number, max_unique_bars)
    end

    function get_rebar_layout_shear_optim(b, zones)
        return 0.
    end

    function get_h(d, zone_layers)

        # Get maximum number of layers and maximum bar size per layer across all zones
        if isempty(zone_layers)
            return d + Av_default + clear_bottom
        end

        non_empty_zones = filter(zone -> !isempty(zone), zone_layers)
        
        if isempty(non_empty_zones)
            return d + Av_default + clear_bottom
        end
        
        # Find zones with maximum number of layers
        n_layers = maximum(length(zone) for zone in non_empty_zones)
        max_layer_zones = filter(zone -> length(zone) == n_layers, non_empty_zones)

        if isempty(max_layer_zones)
            return d + Av_default + clear_bottom
        end
        
        # Among those zones, find the one with largest diameter rebar
        max_zone = nothing
        max_rebar = 0
        
        for zone in max_layer_zones
            zone_max_rebar = maximum(maximum(layer) for layer in zone)
            if zone_max_rebar > max_rebar
                max_rebar = zone_max_rebar
                max_zone = zone
            end
        end

        if isnothing(max_zone)
            return d + Av_default + clear_bottom
        end
        
        # Get maximum diameters for each layer from the controlling zone
        max_rebars = [maximum(layer) for layer in max_zone]
        max_diameters = [REBAR_TYPES[max_rebar].diameter for max_rebar in max_rebars]
    
        # Calculate total rebar height:
        total_rebar_height = 0.5 * max_diameters[1] + 1 * (length(max_diameters) - 1) + sum(max_diameters[2:end-1]) + 0.5 * max_diameters[end]
        h = d + Av_default + 0.5 * total_rebar_height + clear_bottom # only one Av_diameter because the other is accounted for in the initialization
        
        return h
    
    end

    function rc_objective(vars)

        b, d = vars
        println(b, "...", d)
        d += Av_default
        
        for zone in zones
            Mu = ϕb * zone.max_Mn
            zone.As, zone.As′ = get_longitudinal_As_differentiable(b, d, d′, Mu, material)
            if zone.As < 0 || zone.As′ < 0
                println("beam is too small")
            end
        end

        b_rebar = b - 2 * (clear_side + Av_default)
        
        current_As = get_rebar_layout_tension_optim(b_rebar, zones) # max 2 layers for tension
        current_As′ = get_rebar_layout_compression_optim(b_rebar, zones)

        if current_As < 0 || current_As′ < 0
            println("rebar doesn't fit")
        end

        fit_penalty = BIG_M * ((current_As < 0 ? 1 : 0) + (current_As′ < 0 ? 1 : 0))

        steel_volume = (current_As + current_As′) * beam_length
        concrete_volume = b * max_h - steel_volume
        ec = steel_volume * ρ_STEEL_KIPIN3 * ECC_STEEL + concrete_volume * ρ_CONCRETE_KIPIN3 * ECC_CONCRETE

        println(b,"...",d,"...",current_As, "...", current_As′,"...", ec)
        
        objective_value = ec + fit_penalty

        println("objective value: ", objective_value)

        return objective_value

    end

    function rc_constraint(vars)
        
        b, d = vars
        d += Av_default

        h = get_h(d, [zone.rebar for zone in zones])
        h_constraint = h - max_h # h <= max_h

        min_ratio = d / b - 2.5 # d/b <= 2.5
        max_ratio = 1 - d / b # d/b >= 1

        As_constraint = -1 * minimum([zone.As for zone in zones]) # should be positive
        As′_constraint = -1 * minimum([zone.As′ for zone in zones]) # should be positive

        compression_steel_yielding = check_compression_steel_yielded(b, d, d′, current_As, current_As′, fy, fc′)
        tension_controlled = check_tension_controlled_section(b, d, d′, current_As, current_As′, fy, fc′) .* 100

        println([h_constraint, min_ratio, max_ratio, As_constraint, As′_constraint, tension_controlled, compression_steel_yielding])

        return [h_constraint, min_ratio, max_ratio, As_constraint, As′_constraint, tension_controlled, compression_steel_yielding]
    end

    return rc_objective, rc_constraint
end

function size_concrete_zones(zones; material::ReinforcedConcrete=ReinforcedConcrete(5,60), ϕb::Float64 = 0.9, ϕv::Float64 = 0.75, d′=1., clear_side=1., clear_bottom=2.5, clear_top=1., plot::Bool=false)

    objective_function, constraint_function = get_rc_objective(zones, material, max_h=24.0, plot=false)

    lb = [10.0, 10.0]
    ub = [25.0, 30.0]
    x0 = [15.0, 15.0]
    
    model = Nonconvex.Model(objective_function)
    Nonconvex.addvar!(model, lb, ub);
    Nonconvex.add_ineq_constraint!(model, x -> constraint_function(x))

    alg = NLoptAlg(:LN_COBYLA)
    options = NLoptOptions()
    r = optimize(model, alg, ub, options = options)

    return r

end
