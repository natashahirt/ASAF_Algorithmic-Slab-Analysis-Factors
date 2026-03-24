function select_rebar_gurobi(required_As::Float64, b::Float64; clear_side::Float64 = 1.0, max_layers::Int = 2, max_bar_size = 11, bars_per_layer = 6, min_bar_number = 2, max_unique_bars = 2)

    if required_As < 1e-3
        return 0.0, Dict()
    end

    bar_nums = sort([key for key in collect(keys(REBAR_TYPES)) if key <= max_bar_size])
    diameters = Dict(k => v.diameter for (k, v) in REBAR_TYPES)
    areas = Dict(k => v.area for (k, v) in REBAR_TYPES)
    spacing = Dict(i => min(1.0, diameters[i]) for i in bar_nums)

    #model = JuMP.Model(Gurobi.Optimizer)
    model = JuMP.direct_model(Gurobi.Optimizer(gurobi_env)) # reuse it

    JuMP.set_optimizer_attribute(model, "OutputFlag", 0)

    # Set variables
    JuMP.@variable(model, x[bar_nums, 1:max_layers] >= 0, Int)  # number of bars of each type in each layer
    JuMP.@variable(model, y[bar_nums, 1:max_layers] >= 0, Int)  # half the number of bars
    JuMP.@variable(model, z[bar_nums, 1:max_layers], Bin)       # optional centerline bar
    JuMP.@variable(model, total_bar_count[1:max_layers] >= 0, Int) # count of total bars per layer
    JuMP.@constraint(model, [l in 1:max_layers], total_bar_count[l] == sum(x[i, l] for i in bar_nums)) # constraint to make sure it updates

    # No more than 6 bars per layer
    JuMP.@constraint(model, sum(total_bar_count[l] for l in 1:max_layers) >= min_bar_number)
    JuMP.@constraint(model, [l in 1:max_layers], total_bar_count[l] <= bars_per_layer)
    
    # No more than 2 unique bars
    JuMP.@variable(model, bar_used[bar_nums], Bin) # Create binary variable to track if a bar size is used
    M = 100  # Large enough to cover max possible y[i,l]
    for i in bar_nums, l in 1:max_layers
        JuMP.@constraint(model, y[i,l] <= M * bar_used[i])  # if bar_used[i] == 0, y must be 0
        JuMP.@constraint(model, z[i,l] <= bar_used[i])      # if bar is used as centerline, mark as used
    end
    JuMP.@constraint(model, sum(bar_used[i] for i in bar_nums) <= max_unique_bars)

    # Limit unique bars to be no more than 2 standard sizes apart
    JuMP.@variable(model, min_bar_idx >= 1, Int)
    JuMP.@variable(model, max_bar_idx <= length(bar_nums), Int)
    for (j, bar_num) in enumerate(bar_nums)
        JuMP.@constraint(model, min_bar_idx <= j + (1 - bar_used[bar_num]) * length(bar_nums)) # If bar_num is used, then min_bar_idx <= j
        JuMP.@constraint(model, max_bar_idx >= j - (1 - bar_used[bar_num]) * length(bar_nums)) # If bar_num is used, then max_bar_idx >= j
    end
    JuMP.@constraint(model, max_bar_idx - min_bar_idx <= 2)

    # Symmetry constraint
    JuMP.@constraint(model, [i in bar_nums, l in 1:max_layers], x[i, l] == 2 * y[i, l] + z[i, l]) # bars are either symmetric or with a centerline
    JuMP.@constraint(model, [l in 1:max_layers], sum(z[i, l] for i in bar_nums) <= 1) # only one bar type can be symmetric/odd numbered

    # Total As constraint, maximum bars per layer constraint
    tol = 0.1
    JuMP.@constraint(model, sum(areas[i] * x[i, l] for i in bar_nums for l in 1:max_layers) >= required_As - tol)
    JuMP.@constraint(model, sum(areas[i] * x[i, l] for i in bar_nums for l in 1:max_layers) <= required_As + tol)
    for l in 1:max_layers
        JuMP.@constraint(model, sum(x[i, l] for i in bar_nums) <= bars_per_layer) # Limit number of bars per layer to 6
        JuMP.@constraint(model, sum(x[i, l] * (diameters[i] + spacing[i]) for i in bar_nums) - minimum(values(spacing)) <= b - 2 * clear_side)
    end

    # Ensure heaviest bars go in the bottom layer
    for i in bar_nums
        for l in 2:max_layers
            JuMP.@constraint(model, x[i, l] <= x[i, 1])
        end
    end

    # Objective: minimize total steel area and penalize additional layers
    JuMP.@variable(model, layer_used[2:max_layers], Bin) # Binary variable for layer usage
    for l in 2:max_layers
        JuMP.@constraint(model, sum(x[i,l] for i in bar_nums) <= 1000 * layer_used[l])
        JuMP.@constraint(model, sum(x[i,l] for i in bar_nums) >= layer_used[l])
    end

    sorted_bars = sort(bar_nums, rev=true)  # largest to smallest

    for l in 2:max_layers
        for (i, bar_size) in enumerate(sorted_bars)
            sum_bottom = sum(x[s, 1] for s in sorted_bars[1:i]) # sum of bars in bottom layer of size >= bar_sizez
            sum_upper = sum(x[s, l] for s in sorted_bars[1:i]) # sum of bars in upper layer of size >= bar_size
            JuMP.@constraint(model, sum_bottom >= sum_upper)
        end
    end

    # JuMP.@objective(model, Min, 
    #     sum(areas[i] * x[i, l] for i in bar_nums for l in 1:max_layers) + 
    #     10.0 * sum(layer_used[l] for l in 2:max_layers)
    # )

    # penalize the number of bars
    λ = 2
    JuMP.@objective(model, Min, 
        sum(areas[i] * x[i, l] for i in bar_nums for l in 1:max_layers) + 
        10.0 * sum(layer_used[l] for l in 2:max_layers) +
        λ * sum(x[i, l] for i in bar_nums for l in 1:max_layers)
    )

    JuMP.optimize!(model)
    
    # Calculate total steel area from optimal solution    
    if JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        final_As = sum(areas[i] * JuMP.value(x[i,l]) for i in bar_nums for l in 1:max_layers)
        solution = Dict()
        for i in bar_nums, l in 1:max_layers
            n = JuMP.value(x[i, l])
            if n > 0.5
                solution[(i, l)] = round(Int, n)
            end
        end
        return final_As, solution
    else
        error("No feasible bar configuration found.")
    end
end

function select_rebar_gurobi(
        zones::Vector{LongitudinalZone}, 
        beam_ids::Vector{Int64},
        b::Float64; 
        clear_side::Float64 = 1.0, 
        max_layers::Int = 2, 
        max_bar_size = 11, 
        bars_per_layer = 6, 
        min_bar_number = 2, 
        max_unique_bars = 2,
        force::Symbol = :compression
)

    # setup
    bar_nums = sort([key for key in collect(keys(REBAR_TYPES)) if key <= max_bar_size])
    diameters = Dict(k => v.diameter for (k, v) in REBAR_TYPES)
    areas = Dict(k => v.area for (k, v) in REBAR_TYPES)
    spacing = Dict(i => min(1.0, diameters[i]) for i in bar_nums)

    model = JuMP.direct_model(Gurobi.Optimizer(gurobi_env))
    JuMP.set_optimizer_attribute(model, "OutputFlag", 0)

    M = 100  # Large enough to cover max possible y[i,l]

    # get longitudinal zones as well as respective beam ids
    Z = length(zones)

    groups = unique(beam_ids)

    # group vars
    JuMP.@variable(model, global_bar_used[g in groups, i=bar_nums], Bin)
    JuMP.@variable(model, global_min_bar_idx[g in groups], lower_bound=1, Int)
    JuMP.@variable(model, global_max_bar_idx[g in groups], upper_bound=length(bar_nums), Int)

    # individual variables
    JuMP.@variable(model, x[zidx=1:Z, i=bar_nums, l=1:max_layers] >= 0, Int)  # number of bars of each type in each layer in each zone
    JuMP.@variable(model, y[zidx=1:Z, i=bar_nums, l=1:max_layers] >= 0, Int)  # half the number of bars
    JuMP.@variable(model, z[zidx=1:Z, i=bar_nums, l=1:max_layers], Bin)       # optional centerline bar
    JuMP.@variable(model, total_bar_count[zidx=1:Z, l=1:max_layers] >= 0, Int) # count of total bars per layer
    JuMP.@variable(model, bar_used[zidx=1:Z, i=bar_nums], Bin) # Create binary variable to track if a bar size is used per zone
    JuMP.@variable(model, layer_used[zidx=1:Z, l=2:max_layers], Bin) # Binary variable for layer usage

    # individual zone limits
    for (zidx, zone) in enumerate(zones)

        g = beam_ids[zidx] # which group do you belong in?

        # get zone As        
        if force == :tension
            zone_As = zone.As′
        elseif force == :compression
            zone_As = zone.As
        else
            error("force must be :tension or :compression")
        end

        if zone_As < 1e-3
            for l in 1:max_layers
                JuMP.@constraint(model, total_bar_count[zidx, l] == 0)
                for i in bar_nums
                    JuMP.@constraint(model, x[zidx, i, l] == 0)
                    JuMP.@constraint(model, y[zidx, i, l] == 0)
                    JuMP.@constraint(model, z[zidx, i, l] == 0)
                end
            end
            zones[zidx].rebar_top = Vector{Vector{Int64}}()
            zones[zidx].rebar_bottom = Vector{Vector{Int64}}()
            continue
        end

        JuMP.@constraint(model, [l in 1:max_layers], total_bar_count[zidx, l] == sum(x[zidx, i, l] for i in bar_nums)) # constraint to make sure it updates 
        
        # min total bar count across layers
        JuMP.@constraint(model, sum(total_bar_count[zidx, l] for l in 1:max_layers) >= min_bar_number)

        # count unique bars
        for i in bar_nums, l in 1:max_layers
            JuMP.@constraint(model, y[zidx,i,l] <= M * bar_used[zidx,i])  # if bar_used == 0, y must be 0
            JuMP.@constraint(model, z[zidx,i,l] <= bar_used[zidx,i])      # if bar is used as centerline, mark as used
        end
        
        # limit number of unique bars
        JuMP.@constraint(model, sum(bar_used[zidx,i] for i in bar_nums) <= max_unique_bars) 

        # symmetry
        JuMP.@constraint(model, [i in bar_nums, l in 1:max_layers], x[zidx, i, l] == 2 * y[zidx, i, l] + z[zidx, i, l]) # bars are either symmetric or with a centerline
        JuMP.@constraint(model, [l in 1:max_layers], sum(z[zidx, i, l] for i in bar_nums) <= 1) # only one bar type can be symmetric/odd numbered

        # As target
        tol = 0.1
        JuMP.@constraint(model, sum(areas[i] * x[zidx,i,l] for i in bar_nums for l in 1:max_layers) >= zone_As - tol)
        JuMP.@constraint(model, sum(areas[i] * x[zidx,i,l] for i in bar_nums for l in 1:max_layers) <= zone_As + tol)

        for l in 1:max_layers
            JuMP.@constraint(model, sum(x[zidx,i,l] for i in bar_nums) <= bars_per_layer) # Limit number of bars per layer to 6
            JuMP.@constraint(model, sum(x[zidx,i,l] * (diameters[i] + spacing[i]) for i in bar_nums) - minimum(values(spacing)) <= b - 2 * clear_side)
        end

        # ensure heaviest bars go in the bottom layer
        for i in bar_nums
            for l in 2:max_layers
                JuMP.@constraint(model, x[zidx, i, l] <= x[zidx, i, 1])
            end
        end

        # layer usage coupling
        for l in 2:max_layers
            JuMP.@constraint(model, sum(x[zidx, i,l] for i in bar_nums) <= 1000 * layer_used[zidx, l])
            JuMP.@constraint(model, sum(x[zidx, i,l] for i in bar_nums) >= layer_used[zidx, l])
        end

        sorted_bars = sort(bar_nums, rev=true)  # largest to smallest

        for l in 2:max_layers
            for (i, bar_size) in enumerate(sorted_bars)
                sum_bottom = sum(x[zidx, s, 1] for s in sorted_bars[1:i])
                sum_upper = sum(x[zidx, s, l] for s in sorted_bars[1:i])
                JuMP.@constraint(model, sum_bottom >= sum_upper)
            end
        end

        # global bar usage
        for i in bar_nums, l in 1:max_layers
            JuMP.@constraint(model, y[zidx, i, l] <= M * global_bar_used[g, i])
            JuMP.@constraint(model, z[zidx, i, l] <= global_bar_used[g, i])
        end

    end

    # group limits
    for g in groups
        # limit the number of unique bars used
        JuMP.@constraint(model, sum(global_bar_used[g, i] for i in bar_nums) <= max_unique_bars)
        for (j, bar_num) in enumerate(bar_nums)
            JuMP.@constraint(model, global_min_bar_idx[g] <= j + (1 - global_bar_used[g, bar_num]) * length(bar_nums))
            JuMP.@constraint(model, global_max_bar_idx[g] >= j - (1 - global_bar_used[g, bar_num]) * length(bar_nums))
        end
        JuMP.@constraint(model, global_max_bar_idx[g] - global_min_bar_idx[g] <= 2)
    end

    # penalize the number of bars
    λ = 2
    JuMP.@objective(model, Min, 
        sum(areas[i] * x[zidx, i, l] for zidx in 1:Z for i in bar_nums for l in 1:max_layers) + 
        10.0 * sum(layer_used[zidx, l] for zidx in 1:Z for l in 2:max_layers) +
        λ * sum(x[zidx, i, l] for zidx in 1:Z for i in bar_nums for l in 1:max_layers)
    )

    JuMP.optimize!(model)
    
    # calculate total steel area from optimal solution    
    if JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        per_zone_As = [sum(areas[i] * JuMP.value(x[zidx, i, l]) for i in bar_nums for l in 1:max_layers) for zidx in 1:Z]
        total_As = sum(per_zone_As)
        for (zidx, zone) in enumerate(zones)
            # build layers as bar sizes and drop empty layers
            layers = [Int[] for _ in 1:max_layers]
            for l in 1:max_layers
                for i in bar_nums
                    n = Int(round(JuMP.value(x[zidx, i, l])))
                    if n > 0
                        append!(layers[l], fill(i, n))
                    end
                end
            end
            filter!(!isempty, layers)
            zone_rebar = layers
            provided_As = per_zone_As[zidx]

            tension_on_top = zone.sign > 0
            
            if force == :tension
                zone.actual_As = provided_As
                if tension_on_top
                    zone.rebar_top = zone_rebar
                else
                    zone.rebar_bottom = zone_rebar
                end
            elseif force == :compression
                zone.actual_As′ = provided_As
                if tension_on_top
                    zone.rebar_top = zone_rebar
                else
                    zone.rebar_bottom = zone_rebar
                end
            end
        end

        return total_As, per_zone_As
    else
        error("No feasible bar configuration found.")
    end
end

function get_rebar_layout(b::Float64, zones::Vector{LongitudinalZone}, max_layers, compression, d′, ϕb, clear_side, max_bar_size, bars_per_layer, min_bar_number, max_unique_bars)
    
    bar_nums = sort([key for key in collect(keys(REBAR_TYPES)) if key <= max_bar_size])
    # Create dictionaries
    diameters = Dict{Int, Float64}()
    areas = Dict{Int, Float64}()
    spacing = Dict{Int, Float64}()
    
    for (k, v) in REBAR_TYPES
        if k <= max_bar_size
            diameters[k] = v.diameter
            areas[k] = v.area
            spacing[k] = min(1.0, v.diameter)
        end
    end
    
    # Initialize empty model
    jump_model = JuMP.direct_model(Gurobi.Optimizer(gurobi_env))
    JuMP.set_optimizer_attribute(jump_model, "OutputFlag", 0)

    # Global variables for unique bar control
    JuMP.@variable(jump_model, global_bar_used[i=bar_nums], Bin)
    JuMP.@variable(jump_model, global_min_bar_idx, lower_bound=1, Int)
    JuMP.@variable(jump_model, global_max_bar_idx, upper_bound=length(bar_nums), Int)

    JuMP.@variable(jump_model, x[zone=1:length(zones), i=bar_nums, l=1:max_layers], Int)
    JuMP.@variable(jump_model, y[zone=1:length(zones), i=bar_nums, l=1:max_layers], Int)
    JuMP.@variable(jump_model, z[zone=1:length(zones), i=bar_nums, l=1:max_layers], Bin)

    JuMP.@variable(jump_model, total_bar_count[zone=1:length(zones), l=1:max_layers], Int)
    JuMP.@variable(jump_model, layer_used[zone=1:length(zones), l=2:max_layers], Bin)

    # Add non-negativity constraints
    JuMP.@constraint(jump_model, [zone=1:length(zones), i in bar_nums, l=1:max_layers], x[zone, i, l] >= 0)
    JuMP.@constraint(jump_model, [zone=1:length(zones), i in bar_nums, l=1:max_layers], y[zone, i, l] >= 0)
    JuMP.@constraint(jump_model, [zone=1:length(zones), l=1:max_layers], total_bar_count[zone, l] >= 0)

    for (zone_idx, zone) in enumerate(zones)
        zone_As = compression ? zone.As′ : zone.As
        # Force all x, y, z variables to zero for this zone
        if zone_As < 1e-3
            continue
        end

        # Zone-specific constraints
        JuMP.@constraint(jump_model, [l in 1:max_layers], total_bar_count[zone_idx, l] == sum(x[zone_idx, i, l] for i in bar_nums))
        JuMP.@constraint(jump_model, sum(total_bar_count[zone_idx, l] for l in 1:max_layers) >= min_bar_number)
        
        # Link zone variables to global bar usage
        M = 100  # Large enough to cover max possible y[i,l]
        for i in bar_nums, l in 1:max_layers
            JuMP.@constraint(jump_model, y[zone_idx, i, l] <= M * global_bar_used[i])  # if global_bar_used[i] == 0, y must be 0
            JuMP.@constraint(jump_model, z[zone_idx, i, l] <= global_bar_used[i])      # if bar is used as centerline, mark as used
        end

        # Symmetry constraint
        JuMP.@constraint(jump_model, [i in bar_nums, l in 1:max_layers], x[zone_idx, i, l] == 2 * y[zone_idx, i, l] + z[zone_idx, i, l])
        JuMP.@constraint(jump_model, [l in 1:max_layers], sum(z[zone_idx, i, l] for i in bar_nums) <= 1)

        # Total As constraint
        tol = 0.1
        JuMP.@constraint(jump_model, sum(areas[i] * x[zone_idx, i, l] for i in bar_nums for l in 1:max_layers) >= zone_As - tol)
        JuMP.@constraint(jump_model, sum(areas[i] * x[zone_idx, i, l] for i in bar_nums for l in 1:max_layers) <= zone_As + tol)
        
        # Maximum bars per layer constraint
        for l in 1:max_layers
            JuMP.@constraint(jump_model, sum(x[zone_idx, i, l] for i in bar_nums) <= bars_per_layer)
            JuMP.@constraint(jump_model, sum(x[zone_idx, i, l] * (diameters[i] + spacing[i]) for i in bar_nums) - minimum(values(spacing)) <= b - 2 * clear_side)
        end

        # Ensure heaviest bars go in the bottom layer
        for i in bar_nums
            for l in 2:max_layers
                JuMP.@constraint(jump_model, x[zone_idx, i, l] <= x[zone_idx, i, 1])
            end
        end

        # Layer usage constraints
        for l in 2:max_layers
            JuMP.@constraint(jump_model, sum(x[zone_idx, i, l] for i in bar_nums) <= 1000 * layer_used[zone_idx, l])
            JuMP.@constraint(jump_model, sum(x[zone_idx, i, l] for i in bar_nums) >= layer_used[zone_idx, l])
        end

        # Bar ordering constraints
        sorted_bars = sort(bar_nums, rev=true)
        for l in 2:max_layers
            for (i, bar_size) in enumerate(sorted_bars)
                sum_bottom = sum(x[zone_idx, s, 1] for s in sorted_bars[1:i])
                sum_upper = sum(x[zone_idx, s, l] for s in sorted_bars[1:i])
                JuMP.@constraint(jump_model, sum_bottom >= sum_upper)
            end
        end
    end

    # Global constraints on unique bars
    JuMP.@constraint(jump_model, sum(global_bar_used[i] for i in bar_nums) <= max_unique_bars)

    # Limit unique bars to be no more than 2 standard sizes apart globally
    for (j, bar_num) in enumerate(bar_nums)
        JuMP.@constraint(jump_model, global_min_bar_idx <= j + (1 - global_bar_used[bar_num]) * length(bar_nums))
        JuMP.@constraint(jump_model, global_max_bar_idx >= j - (1 - global_bar_used[bar_num]) * length(bar_nums))
    end
    JuMP.@constraint(jump_model, global_max_bar_idx - global_min_bar_idx <= 2)

    # Global objective: minimize total steel area across all zones and penalize additional layers
    λ = 2
    JuMP.@objective(jump_model, Min, 
        sum(areas[i] * x[zone_idx, i, l] for zone_idx in 1:length(zones) for i in bar_nums for l in 1:max_layers) + 
        10.0 * sum(layer_used[zone_idx, l] for zone_idx in 1:length(zones) for l in 2:max_layers) +
        λ * sum(x[zone_idx, i, l] for zone_idx in 1:length(zones) for i in bar_nums for l in 1:max_layers)
    )

    JuMP.optimize!(jump_model)
    
    if JuMP.termination_status(jump_model) == JuMP.MOI.OPTIMAL
        # Return array of steel areas for each zone
        total_As = 0
        for (zone_idx, zone) in enumerate(zones)
            zone_As = sum(areas[i] * JuMP.value(x[zone_idx, i, l]) for i in bar_nums for l in 1:max_layers)
            total_As += zone_As
        end

        """# Get rebar configuration for each zone
        zone_layers = []
        for zone_idx in 1:length(zones)
            n_layers = max_layers
            layer_sizes = [Int[] for _ in 1:n_layers]
            
            # For each bar size and layer, add bars to layer_sizes
            for i in bar_nums, l in 1:max_layers
                num_bars = Int(JuMP.value(x[zone_idx, i, l]))
                if num_bars > 0
                    append!(layer_sizes[l], fill(i, num_bars))
                end
            end
            
            # Sort and arrange each layer's bars
            for i in 1:lastindex(layer_sizes)
                if !isempty(layer_sizes[i])
                    sorted = sort(layer_sizes[i], rev=true)
                    n = length(sorted)
                    reordered = similar(sorted)
                    
                    if n == 1
                        reordered[1] = sorted[1]
                    else
                        # Find unique bar for middle if odd number
                        unique_bar = nothing
                        if n % 2 == 1
                            for bar in sorted
                                if count(==(bar), sorted) == 1
                                    unique_bar = bar
                                    break
                                end
                            end
                        end
                        
                        remaining = filter(x -> x != unique_bar, sorted)
                        mid_idx = div(n+1, 2)
                        
                        # Fill outer positions symmetrically
                        for (j, bar) in enumerate(remaining)
                            if j % 2 == 1
                                reordered[div(j+1,2)] = bar
                            else
                                reordered[n-div(j,2)+1] = bar
                            end
                        end
                        
                        # Place unique bar in middle if exists
                        if !isnothing(unique_bar)
                            reordered[mid_idx] = unique_bar
                        end
                    end
                    layer_sizes[i] = reordered
                end
            end
            
            # Remove empty layers
            filter!(!isempty, layer_sizes)
            push!(zone_layers, layer_sizes)
        end"""
        
        return total_As
    else
        # Return zeros if optimization fails
        return -1 # [[] for zone in zones]
    end
end

Zygote.@nograd get_rebar_layout

function get_rebar_layers(rebar_config::Dict)
    if isempty(rebar_config)
        return []
    end

    n_layers = maximum(last(bar_layer) for (bar_layer, _) in rebar_config)
    layer_sizes = [Int[] for _ in 1:n_layers]

    # For each bar configuration, add the bar size num_bars times to that layer
    for ((size, layer), num_bars) in rebar_config
        append!(layer_sizes[layer], fill(size, num_bars))
    end
    # Sort and arrange each layer's bars
    for i in 1:lastindex(layer_sizes)
        if !isempty(layer_sizes[i])
            sorted = sort(layer_sizes[i], rev=true)
            n = length(sorted)
            reordered = similar(sorted)
            
            if n == 1
                # Single bars go in the middle
                reordered[1] = sorted[1]
            else
                # For odd number of bars, find the unique bar (if any) to place in middle
                unique_bar = nothing
                if n % 2 == 1
                    for bar in sorted
                        if count(==(bar), sorted) == 1
                            unique_bar = bar
                            break
                        end
                    end
                end
                
                # Place bars symmetrically from outside in
                remaining = filter(x -> x != unique_bar, sorted)
                mid_idx = div(n+1, 2)
                
                # Fill outer positions symmetrically
                for (j, bar) in enumerate(remaining)
                    if j % 2 == 1
                        reordered[div(j+1,2)] = bar
                    else
                        reordered[n-div(j,2)+1] = bar
                    end
                end
                
                # Place unique bar in middle if it exists
                if unique_bar !== nothing
                    reordered[mid_idx] = unique_bar
                end
            end
            layer_sizes[i] = reordered
        end
    end

    return layer_sizes
end

