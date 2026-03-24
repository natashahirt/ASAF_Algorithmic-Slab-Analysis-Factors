"""
shear reinforcement design
"""

function get_Vc(fc′, b, d; λ = 1.0)
    # λ = lightweight reduction factor, 1.0 for normal concrete
    # get the shear resistance of the concrete material
    Vc = λ * 2 * safe_sqrt(fc′ * 1e3) * b * d / 1e3
    return Vc
end

function get_Av(Vu, b, d, fc′, fy; ϕ = 0.75, α = pi/2)
    Vc = ϕ * get_Vc(fc′, b, d)
    fc′_psi = fc′ * 1e3
    fy = min(fy, 60) # ACI code limits yield strength of reinforcement to 60ksi
    sin_cos = sin(α) + cos(α)
    if Vu <= (ϕ * Vc / 2) # low shear demand, no stirrups (also no stirrups in special configuration beams)
        Av = 0
    elseif Vu <= (ϕ * Vc) # moderate shear demand, stirrups for crack control not strength
        s = α == pi/2 ? minimum([d/2, 24]) : minimum([3*d/4, 24])
        # Av = 50 * b * s / (fy * 1e3) # ACI minimum stirrups (50 lb/in²)
        Av = 0.75 * safe_sqrt(fc′_psi) * b * s / (fy * 1e3) # ACI minimum stirrups (Eq 11-13)
    else # high shear demand, stirrups for strength
        # ACI shear limits use sqrt(f'c) with f'c in psi; model inputs store f'c in ksi.
        if (Vu - ϕ * Vc) <= ϕ * 4 * safe_sqrt(fc′_psi) * b * d / 1e3
            if α == pi/2
                s = minimum([d/2, 24])
            else
                s = minimum([3*d/4, 24])
            end
        else
            s = minimum([d/4, 12])
        end
        s = maximum([s, 4])
        Av = (Vu - ϕ * Vc) * s / (ϕ * fy * d * sin_cos)
    end

    return Av

end

function get_stirrup_spacing(Vu, Av, b, d, fc′, fy; ϕ=0.75, α=π/2)
    Vc = ϕ * get_Vc(fc′, b, d)
    fc′_psi = fc′ * 1e3
    sin_cos = sin(α) + cos(α)
    if abs(Vu) <= (ϕ * Vc / 2)
        return 0
    elseif abs(Vu) <= (ϕ * Vc)
        return α == pi/2 ? minimum([d/2, 24]) : minimum([3*d/4, 24])
    else
        if (abs(Vu) - ϕ * Vc) <= ϕ * 4 * safe_sqrt(fc′_psi) * b * d / 1e3
            if α == pi/2
                s_max = minimum([Av * fy * 1e3 / (0.75 * safe_sqrt(fc′_psi) * b), Av * fy / (50 * b), d/2, 24])
            else
                s_max = minimum([3*d/4, 24])
            end
        else
            s_max = minimum([d/4, 12])
        end
        s_max = maximum([s_max, 4]) # spacing less than 4" is undesirable
        s_req = (Av * ϕ * fy * d * sin_cos) / (abs(Vu) - ϕ * Vc)
        if α != pi/2
            s_horiz_max = (3/8) * d * cot(α)
            s_horiz_limit = s_horiz_max / cos(α)
            return min(s_req, s_max, s_horiz_limit)
        end
        return min(s_req, s_max)
    end
end

function get_stirrup_zones(spacing::Vector{Real})
    zones = Vector[]
    for (i,s) in enumerate(spacing)
        if !isempty(zones) && zones[end][1] == s
            zones[end][2][2] = i
        else
            push!(zones, [s, [i,i]])
        end
    end
    return zones
end

function get_stirrup_positions(spacing, x, h, d′, clear_bottom, Vu; α=π/2)
    zones = StirrupZone[]
    all_stirrups = Float64[]

    beam_start = x[1]
    beam_end = x[end]
    stirrup_height = h - d′ - clear_bottom
    dx = α == π/2 ? 0. : stirrup_height / tan(α)
    last_spacing = 0.

    i = 1

    while i <= lastindex(spacing)
        s = spacing[i]
        start_idx = i
    
        while i < length(spacing) && abs(s - spacing[i+1]) < 1e-3
            i += 1
        end
        end_idx = i
    
        x_start = x[start_idx]
        x_end = x[end_idx]
    
        if s < 1e-3
            push!(zones, StirrupZone(s, start_idx:end_idx, (x_start, x_end), Float64[]))
            i += 1
            continue
        end
    
        usable_length = x_end - x_start
        n_spaces = floor(Int, usable_length / s)
        n_stirrups = n_spaces + 1
        leftover = usable_length - n_spaces * s
        offset = leftover / 2
        x1 = x_start + offset
        zone_stirrups = [x1 + (j-1)*s for j in 1:n_stirrups if x1 + (j-1)*s <= x_end + 1e-6]
    
        # Enforce minimum spacing from previous zone
        if !isempty(zones) && !isempty(zone_stirrups)
            idx = findmin(abs.(x .- zone_stirrups[1]))[2]
            shear_at_boundary = Vu[idx]
    
            min_spacing = shear_at_boundary < 0 ? min(last_spacing, s) : max(last_spacing, s)
    
            gap = !isempty(zones[end].stirrups) ? zone_stirrups[1] - zones[end].stirrups[end] : Inf
            if gap < min_spacing
                shift = min_spacing - gap
                zone_stirrups = [x + shift for x in zone_stirrups if x + shift <= x_end + 1e-6]
            end
        end
    
        # Filter stirrups by beam bounds, considering inclination
        stirrups = Float64[]
        for x_pos in zone_stirrups
            sign_shear = 1
            if Vu !== nothing
                idx = findmin(abs.(x .- x_pos))[2]
                sign_shear = sign(Vu[idx])
            end
            dx = α == π/2 ? 0.0 : stirrup_height / tan(α)
            x1_local = min(x_pos, x_pos + dx * (sign_shear <= 0 ? 1 : -1))
            x2_local = max(x_pos, x_pos + dx * (sign_shear <= 0 ? 1 : -1))
    
            if x1_local ≥ beam_start - 1e-6 && x2_local ≤ beam_end + 1e-6
                push!(stirrups, x_pos)
            end
        end
    
        push!(zones, StirrupZone(s, start_idx:end_idx, (x_start, x_end), stirrups))
        append!(all_stirrups, stirrups)
    
        if !isempty(stirrups)
            last_spacing = s
        end
    
        i += 1
    end

    # Globally center all stirrups
    if !isempty(all_stirrups)
        endpoints = []

        for x_pos in all_stirrups
            sign_shear = 1
            if Vu !== nothing
                idx = findmin(abs.(x .- x_pos))[2]
                sign_shear = sign(Vu[idx])
            end
            dx = α == π/2 ? 0.0 : stirrup_height / tan(α)
            x_end = sign_shear <= 0 ? x_pos + dx : x_pos - dx
            push!(endpoints, (x_pos, x_end))
        end

        all_min = minimum(min(x1, x2) for (x1, x2) in endpoints)
        all_max = maximum(max(x1, x2) for (x1, x2) in endpoints)
        stirrup_span = all_max - all_min
        beam_span = x[end] - x[1]

        offset = x[1] + (beam_span - stirrup_span)/2 - all_min
        centered_stirrups = [x_pos + offset for x_pos in all_stirrups]
        
        # Update each zone's stirrups with the centered positions
        stirrup_idx = 1
        for zone in zones
            n_zone_stirrups = length(zone.stirrups)
            if n_zone_stirrups > 0
                zone.stirrups = centered_stirrups[stirrup_idx:stirrup_idx+n_zone_stirrups-1]
                stirrup_idx += n_zone_stirrups
            end
        end
    end

    """for zone in zones
        println("Spacing: ", zone.spacing)
        println("From x = ", zone.x_range[1], " to ", zone.x_range[2])
        println("Stirrups at: ", zone.stirrups)
    end"""

    return zones
end

function print_stirrup_zones(stirrup_zones)
    println("Stirrup Zones:")
    println("-------------------------------------------------------------")
    println("Zone | Spacing | Index Range | X Range         | # Stirrups")
    println("-------------------------------------------------------------")
    for (i, zone) in enumerate(stirrup_zones)
        spacing = zone.spacing
        idx_range = zone.index_range
        x_range = zone.x_range
        n_stirrups = length(zone.stirrups)
        println("$(lpad(i,4)) | $(lpad(round(spacing,digits=2),7)) | [$(lpad(idx_range[1],3)), $(lpad(idx_range[end],3))]   | [$(lpad(round(x_range[1],digits=2),8)), $(lpad(round(x_range[end],digits=2),8))] | $(lpad(n_stirrups,10))")
            i, spacing, idx_range[1], idx_range[2], x_range[1], x_range[2], n_stirrups
    end
    println("-------------------------------------------------------------")
end

function get_stirrup_angle(Vx, Vy)
    θ_crack = atan(Vy, Vx)           # crack direction (principal tension)
    θ_stirrup = θ_crack + π/2        # stirrups go perpendicular
    return θ_stirrup                 # radians
end

function distribute_stirrups(zones, x)
    for zone in zones
        distance_interval = [x[zone[2][1]], x[zone[2][2]]]
        distance = distance_interval[2] - distance_interval[1]
        if zone[1] < 1e-3
            n_stirrups = 0
        else 
            n_stirrups = Int64(ceil((distance) / zone[1]))
        end
        push!(zone, distance_interval)
        push!(zone, n_stirrups)
    end
    return zones
end
