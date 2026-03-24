function get_embodied_carbon_rc(b, h, l, clear_top, clear_bottom, clear_side, As, As′, Av, n_stirrups, α)
    total_beam_volume = b * h * l
    longitudinal_steel_volume = (As + As′) * l
    # Calculate stirrup length based on beam dimensions and cover
    inner_rectangle_height = (h - clear_top - clear_bottom)
    inner_rectangle_width = b - 2*clear_side
    rectangle_length =  (2 * (inner_rectangle_height + inner_rectangle_width)) * 1.05
    stirrup_steel_volume = Av * rectangle_length
    steel_volume = longitudinal_steel_volume + n_stirrups * stirrup_steel_volume
    concrete_volume = total_beam_volume - steel_volume
    concrete_ec = ECC_CONCRETE * ρ_CONCRETE_KIPIN3 * concrete_volume
    steel_ec = ECC_STEEL * ρ_STEEL_KIPIN3 * steel_volume
    ec = concrete_ec + steel_ec
    # Convert from kip to kg
    ec = ec * 1/kg_to_kip(1)
    return ec
end

function get_c(b, d, d′, As, As′, fy, fc′, β1)
    # Calculate neutral axis depth c by solving equilibrium equation
    # Cc + Cs = T
    # 0.85*fc′*β1*b*c + As′*fy*(c-d′)/c = As*fy
    
    a = (As - As′) * fy / (0.85 * fc′ * b)
    c = a / β1
    return c
end

function calculate_height(d::Real, clear_bottom, Av_diameter, layer_sizes)
    maximum_diameters = [maximum([REBAR_TYPES[bar].diameter for bar in layer]) for layer in layer_sizes] # tensile reinforcement
    total_rebar_height = 0.5 * maximum_diameters[1] + 1 * (length(maximum_diameters) - 1) + sum(maximum_diameters[2:end-1]) + 0.5 * maximum_diameters[end]
    h = d + clear_bottom + total_rebar_height / 2 + 2 * Av_diameter
    return h
end
