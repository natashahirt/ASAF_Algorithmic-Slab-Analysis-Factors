# everything that goes into this function should be in inches and kips
function size_concrete_section(b, d, x::Vector{Float64}, Mn::Vector{Float64}, Vn::Vector{Float64}, fc′::Float64, fy::Float64; ϕb::Float64 = 0.9, ϕv::Float64 = 0.75, d′=1., clear_side=1., clear_bottom=2.5, clear_top=1., α=π/2, use_critical_Vn::Bool=true, rounded=false, plot::Bool=false, verbose::Bool=true)
    #initialize outputs
    final_As = final_As′ = final_Av = 0
    
    # get the design variables
    max_beam_Mn = maximum(abs.(Mn))
    Mu = ϕb * max_beam_Mn

    # in case we want to use critical Vn instead of full Vn
    if use_critical_Vn
        Vn = copy(Vn)
        left_critical_idx = findfirst(x -> x >= d, x)
        right_critical_idx = findlast(x -> x <= (maximum(x) - d), x)
        if !isnothing(left_critical_idx) && !isnothing(right_critical_idx)
            fill!(@view(Vn[1:left_critical_idx]), Vn[left_critical_idx])
            fill!(@view(Vn[right_critical_idx:end]), Vn[right_critical_idx])
        end
    end

    max_beam_Vn = maximum(abs.(Vn))
    Vu = ϕv * max_beam_Vn
    Vu_vec = ϕv * Vn
    
    # find the stirrup size
    Av = get_Av(Vu, b, d, fc′, fy, α=α)
    required_Av, rebar_config_V = select_rebar_gurobi(Av, b, bars_per_layer=4, max_layers=1, min_bar_number=2)
    layer_sizes_V = get_rebar_layers(rebar_config_V)[1] # since it should only be one layer remove the outside list
    Av_diameter, Av_area = REBAR_TYPES[maximum(layer_sizes_V)].diameter, REBAR_TYPES[maximum(layer_sizes_V)].area * length(layer_sizes_V)
    final_Av = length(layer_sizes_V)*Av_diameter

    nom_clear_top = required_As = required_As′ = 0.
    layer_sizes = layer_sizes′ = [[]]

    while nom_clear_top < clear_top && d′ <= 3
        # find the tensile reinforcement
        ρ_single_required = singly_reinforced_required_ρ(Mu, b, d, fy, fc′)
        ρ_single_max = singly_reinforced_max_ρ(fc′, fy)
        ρ_single_min = singly_reinforced_min_ρ(fc′, fy)
        ρ_estimated = (Mu/(0.9*d*fy)) / (b * d)

        if verbose
            println("ρ max: ", ρ_single_max)
            println("ρ min: ", ρ_single_min)
            println("ρ required: ", ρ_single_required)
            println("ρ estimated: ", ρ_estimated)
        end

        # find the compression reinforcement
        if ρ_single_max < ρ_single_required
            Mu₂ = get_Mu₂(ρ_single_max, b, d, fy, fc′)
            verbose && println("Mu₂: ", Mu₂)
            Mu₁, required_As′ = get_Mu₁(Mu, Mu₂, d, d′, fy)
            verbose && println("Mu₁: ", Mu₁)
            required_As′ = maximum([required_As′, 0.0])
            required_As = ρ_single_max * b * d + required_As′
        elseif ρ_single_min > ρ_single_required
            required_As′ = 0.0
            required_As = ρ_single_min * b * d
        else
            required_As′ = 0.0
            required_As = ρ_single_required * b * d
        end

        if verbose
            println("As required: ", required_As)
            println("As permitted: ", minimum([ρ_single_max, ρ_single_required]) * b * d)
            println("As′ required: ", required_As′)
        end

        # get the arrangement and required steel sizes
        final_As, rebar_config = select_rebar_gurobi(required_As, b-2*(clear_side+Av_diameter))
        final_As′, rebar_config′ = select_rebar_gurobi(required_As′, b-2*(clear_side+Av_diameter), max_layers=1)
        layer_sizes = get_rebar_layers(rebar_config)
        layer_sizes′ = get_rebar_layers(rebar_config′)

        verbose && println("As reinforced: ", final_As)

        if required_As′ > 0
            nom_clear_top = d′ - REBAR_TYPES[maximum(layer_sizes′[1])].diameter / 2 - Av_diameter
        else
            nom_clear_top = d′ - Av_diameter
        end
        verbose && println("clear top: ", nom_clear_top)
        if nom_clear_top < clear_top
            d′ += 0.1
        end
        if verbose
            println("Layer sizes:")
            println(layer_sizes)
            println("Layer sizes′:")
            println(layer_sizes′)
        end
    end

    h = calculate_height(d, clear_bottom, Av_diameter, layer_sizes)
    spacing = [get_stirrup_spacing(V, Av_area, b, d, fc′, fy, α=α) for V in Vu_vec]
    stirrup_positions = get_stirrup_positions(spacing, x, h, d′, clear_bottom, Vu_vec, α=α)

    if plot
        # find the stirrup spacing and analyze zones using the new comprehensive function
        #stirrup_zones = get_stirrup_zones(spacing)
        #stirrup_zones = distribute_stirrups(stirrup_zones, x)
        fig_section = plot_beam_section(b, h, d, clear_bottom, clear_side, d′ + Av_diameter, layer_sizes, layer_sizes′, Av_diameter, rounded=rounded)
        fig_stirrups = plot_stirrups(x, stirrup_positions, h, d′, clear_bottom; α=α, shear_signs=Vn)
        display(fig_section)
        # display(fig_stirrups)
    end

    n_stirrups = sum([length(zone.stirrups) for zone in stirrup_positions])
    EC = get_embodied_carbon_rc(b, h, l, clear_top, clear_bottom, clear_side, final_As, final_As′, Av_diameter, n_stirrups, α)
    verbose && println(EC)

    # checks
    check_compression = check_compression_steel_yielded(b, d, d′, final_As, final_As′, fy, fc′)
    check_tension = check_tension_controlled_section(b, d, d′, final_As, final_As′, fy, fc′)
    if verbose
        println(check_compression <= 0 ? "Compression steel has yielded" : "Compression steel has not yielded")
        println(check_tension <= 0 ? "Tension controlled section" : "Not tension controlled section")
    end

    return final_As, final_As′, final_Av, h, d′
end
