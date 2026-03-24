
function plot_beam_section(b, h, d, clear_bottom, clear_side, d′, layer_sizes, layer_sizes′, Av_diameter; rounded=false)
    # Create figure
    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())
    
    if rounded
        h = ceil(h)
    end

    # Plot concrete outline
    concrete = Rect(-b/2, 0, b, h)
    poly!(ax, concrete, color=:lightgray)
    lines!(ax, [-b/2,b/2], [h-d, h-d], color=:black, linestyle=:dash)

    # Plot rebar
    # Get number of layers by taking maximum layer number from rebar config
    # rebar_config is Dict{Any,Any} mapping (bar_size, layer) => num_bars
    # Get number of layers and initialize arrays for each layer
    centers = compute_bar_centers(layer_sizes, b, h, d, clear_side, Av_diameter)
    min_y = h
    
    for layer in centers # Plot circles for rebar
        for bar in layer
            x, y, diameter = bar
            min_y = minimum([min_y, y - diameter/2])
            circle = Circle(Point2f(x, y), diameter/2)
            poly!(ax, circle, color=:black)
        end
    end

    max_y = h - d′
    if !isempty(layer_sizes′)
        centers′ = compute_bar_centers(layer_sizes′, b, h, d′, clear_side, Av_diameter)
        for layer in centers′ # Plot circles for rebar
            for bar in layer
                x, y, diameter = bar
                max_y = maximum([max_y, y + diameter / 2])
                circle = Circle(Point2f(x, y), diameter/2)
                poly!(ax, circle, color=:black)
            end
        end
        text!(ax, "d′ = $(round(d′,digits=2))\"", position=(b/2+1, h-d′), rotation=pi/2, align=(:center, :center))
        lines!(ax, [-b/2,b/2], [h-d′, h-d′], color=:black, linestyle=:dash)
    end

    # Plot stirrups
    # Get coordinates of inner rectangle corners
    inner_x1 = -b/2 + clear_side + Av_diameter
    inner_x2 = b/2 - clear_side - Av_diameter
    inner_y1 = min_y
    inner_y2 = max_y
    
    # Draw inner rectangle
    lines!(ax, [inner_x1, inner_x2], [inner_y1, inner_y1], color=:black, linestyle=:dot)
    lines!(ax, [inner_x2, inner_x2], [inner_y1, inner_y2], color=:black, linestyle=:dot)
    lines!(ax, [inner_x2, inner_x1], [inner_y2, inner_y2], color=:black, linestyle=:dot)
    lines!(ax, [inner_x1, inner_x1], [inner_y2, inner_y1], color=:black, linestyle=:dot)

    # Get coordinates of outer rectangle corners
    outer_x1 = inner_x1 - Av_diameter
    outer_x2 = inner_x2 + Av_diameter
    outer_y1 = inner_y1 - Av_diameter
    outer_y2 = inner_y2 + Av_diameter
    
    # Draw outer rectangle
    lines!(ax, [outer_x1, outer_x2], [outer_y1, outer_y1], color=:black, linestyle=:dot)
    lines!(ax, [outer_x2, outer_x2], [outer_y1, outer_y2], color=:black, linestyle=:dot)
    lines!(ax, [outer_x2, outer_x1], [outer_y2, outer_y2], color=:black, linestyle=:dot)
    lines!(ax, [outer_x1, outer_x1], [outer_y2, outer_y1], color=:black, linestyle=:dot)
    
    # Add dimensions with more padding for text
    text!(ax, "b = $(round(b,digits=2))\"", position=(0, -1), rotation=0, align=(:center, :center))
    text!(ax, "h = $(round(h,digits=2))\"", position=(b/2+1, h/2), rotation=pi/2, align=(:center, :center))
    text!(ax, "d = $(round(d,digits=2))\"", position=(b/2+1, h-d), rotation=pi/2, align=(:center, :center))

    # Set axis limits with more padding to show text
    xlims!(ax, -b/2-6, b/2+6)
    ylims!(ax, -2, h+2)
    
    hidedecorations!(ax)
    hidespines!(ax)
    
    return fig
end

function compute_bar_centers(layers, b, h, d, clear_side, Av_diameter)
    diameters = [[REBAR_TYPES[bar].diameter for bar in layer] for layer in layers]
    maximum_diameters = [maximum(diameters[i]) for (i, layer) in enumerate(layers)]
    bar_centers = []
    
    # Start from bottom and work up
    total_height = sum(maximum_diameters) + (length(layers)-1) * 1.0 # Total height of all layers including spacing
    start_y = h - d - total_height/2 # Center around h-d
    
    for (i, layer) in enumerate(layers)
        # Calculate y position working up from start_y
        y = start_y + (sum(maximum_diameters[1:i-1]) + (i-1) * 1.0 + 0.5 * maximum_diameters[i]) 
        
        n_bars = length(layer)
        total_bars_width = sum(diameters[i])
        available_width = b - 2 * clear_side - 2 * Av_diameter
        
        if n_bars == 1
            push!(bar_centers, [(0, y, diameters[i][1])])
            continue
        end
        
        total_spacing = available_width - total_bars_width
        spacing = total_spacing / (n_bars - 1)
        
        xs = [clear_side + Av_diameter + diameters[i][1]/2]
        for j in 2:n_bars
            x_new = xs[end] + diameters[i][j-1]/2 + spacing + diameters[i][j]/2 
            push!(xs, x_new)
        end
        
        push!(bar_centers, [(x - b/2, y, d) for (x, d) in zip(xs, diameters[i])])
    end
    
    return bar_centers
end

function plot_stirrups(x_inches, stirrup_zones, h, d′, clear_bottom; α=π/2, shear_signs=nothing)
    fig = Figure();
    ax = Axis(fig[1,1], aspect=DataAspect(), xlabel="Length (in)", ylabel="Height (in)");
    xlims!(ax, [0,x_inches[end]])
    
    # concrete base
    concrete = Rect(0, 0, x_inches[end], h)
    poly!(ax, concrete, color=:lightgray)

    all_stirrups = centered_stirrups_with_min_gap(stirrup_zones, shear_signs, x_inches)
    all_stirrups = center_stirrups_with_angles(all_stirrups, x_inches, h, d′, clear_bottom, α, shear_signs)

    for stirrup in all_stirrups
        # Determine direction based on shear sign
        local_sign = 1
        if shear_signs !== nothing
            # Find the closest x in x_inches to stirrup
            idx = findmin(abs.(x_inches .- stirrup))[2]
            local_sign = sign(shear_signs[idx])
        end
        if α == π/2
            # Vertical stirrup
            lines!(ax, [stirrup, stirrup], [clear_bottom, h-d′], color=:black)
        else
            stirrup_height = h - d′ - clear_bottom
            dx = stirrup_height / tan(α)
            inclined_x = local_sign <= 0 ? stirrup + dx : stirrup - dx
            if min(stirrup, inclined_x) >= 0 && max(stirrup, inclined_x) <= x_inches[end]
                lines!(ax, [stirrup, inclined_x], [clear_bottom, h - d′], color=:black)
            end
        end
    end
    return fig
end

function plot_stirrups(x_inches, zones, h, d′, clear_bottom; α=π/2, shear_signs=nothing)
    fig = Figure();
    ax = Axis(fig[1,1], aspect=DataAspect(), xlabel="Length (in)", ylabel="Height (in)");
    xlims!(ax, [0,x_inches[end]])
    
    # concrete base
    concrete = Rect(0, 0, x_inches[end], h)
    poly!(ax, concrete, color=:lightgray)

    for zone in zones
        for stirrup in zone.stirrups
            local_sign = 1
            if shear_signs !== nothing
                idx = findmin(abs.(x_inches .- stirrup))[2]
                local_sign = sign(shear_signs[idx])
            end
    
            if α == π/2
                # Vertical stirrup
                lines!(ax, [stirrup, stirrup], [clear_bottom, h - d′], color=:black)
            else
                stirrup_height = h - d′ - clear_bottom
                dx = stirrup_height / tan(α)
                x_end = local_sign <= 0 ? stirrup + dx : stirrup - dx
                if min(stirrup, x_end) >= 0 && max(stirrup, x_end) <= x_inches[end]
                    lines!(ax, [stirrup, x_end], [clear_bottom, h - d′], color=:black)
                end
            end
        end
    end

    return fig
end

function centered_stirrups_with_min_gap(zones, shear_signs, x_inches)
    all_stirrups = Float64[]
    last_stirrup = nothing
    last_spacing = nothing

    for (i, zone) in enumerate(zones)
        spacing = zone[1]
        x_start, x_end = zone[3]
        if spacing < 1e-3
            continue
        end
        # Centered positions in this zone
        n_spaces = floor((x_end - x_start) / spacing)
        n_stirrups = Int(n_spaces) + 1
        leftover = (x_end - x_start) - n_spaces * spacing
        offset = leftover / 2
        x1 = x_start + offset
        x_stirrups = [x1 + (j-1)*spacing for j in 1:n_stirrups if x1 + (j-1)*spacing <= x_end + 1e-6]

        # Enforce minimum spacing from previous zone
        if !isempty(all_stirrups) && !isempty(x_stirrups)
            # Find the shear at the boundary (at the first stirrup of this zone)
            idx = findmin(abs.(x_inches .- x_stirrups[1]))[2]
            shear_at_boundary = shear_signs[idx]
            if shear_at_boundary < 0
                min_spacing = min(last_spacing, spacing)
            else
                min_spacing = max(last_spacing, spacing)
            end
            gap = x_stirrups[1] - all_stirrups[end]
            if gap < min_spacing
                shift = min_spacing - gap
                x_stirrups = [x + shift for x in x_stirrups if x + shift <= x_end + 1e-6]
            end
        end

        append!(all_stirrups, x_stirrups)
        if !isempty(x_stirrups)
            last_stirrup = x_stirrups[end]
            last_spacing = spacing
        end
    end
    return all_stirrups
end

function center_stirrups_with_angles(all_stirrups, x_inches, h, d′, clear_bottom, α=π/2, shear_signs=nothing)
    # Compute dx for each stirrup
    stirrup_height = h - d′ - clear_bottom
    endpoints = []
    for (i, x) in enumerate(all_stirrups)
        local_sign = 1
        if shear_signs !== nothing
            idx = findmin(abs.(x_inches .- x))[2]
            local_sign = sign(shear_signs[idx])
        end
        if α == π/2
            push!(endpoints, (x, x))
        else
            dx = stirrup_height / tan(α)
            x_end = local_sign <= 0 ? x + dx : x - dx
            push!(endpoints, (x, x_end))
        end
    end
    # Find min and max of all endpoints
    span_min = minimum([min(x1, x2) for (x1, x2) in endpoints])
    span_max = maximum([max(x1, x2) for (x1, x2) in endpoints])
    stirrup_span = span_max - span_min
    beam_span = x_inches[end] - x_inches[1]
    offset = (beam_span - stirrup_span) / 2 - (span_min - x_inches[1])
    # Shift all stirrups
    all_stirrups_centered = [x + offset for x in all_stirrups]
    return all_stirrups_centered
end