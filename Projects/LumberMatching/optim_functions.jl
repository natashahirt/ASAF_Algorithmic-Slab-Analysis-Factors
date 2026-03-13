function match_segments(cost_matrix::Matrix{Float64}, singleton_vector::Vector{Float64})
    # cost_matrix is symmetric, 1e9 means infeasible
    n = size(cost_matrix, 1)
    model = JuMP.Model(Gurobi.Optimizer)
    JuMP.set_silent(model)

    # Decision variables
    JuMP.@variable(model, x[i=1:n-1, j=i+1:n], Bin)  # x[i,j] = 1 if i and j are matched
    JuMP.@variable(model, z[1:n], Bin)              # z[i] = 1 if i is unmatched (singleton)

    # Each segment appears in exactly one match or singleton
    for i in 1:n
        JuMP.@constraint(model,
            sum(x[min(i,j), max(i,j)] for j in 1:n if i != j && min(i,j) < max(i,j)) + z[i] == 1
        )
    end

    # Objective: minimize excess material
    JuMP.@objective(model, Min,
        sum(cost_matrix[i,j] * x[i,j] for i in 1:n-1, j in i+1:n) +
        sum(singleton_vector[i] * z[i] for i in 1:n)
    )

    JuMP.optimize!(model)

    # Extract results
    matches = []
    singles = []
    for i in 1:n-1, j in i+1:n
        if JuMP.value(x[i,j]) > 0.5
            push!(matches, (i, j))
        end
    end
    
    # Parse singletons and infeasibles
    for i in 1:n
        if JuMP.value(z[i]) > 0.5
            push!(singles, i)
        end
    end

    return matches, singles, JuMP.objective_value(model)
end

function get_costs(segments; board_depths=Float64[])
    singleton_cost, singleton_boards = get_singleton_cost_vector(segments, board_depths=board_depths)

    n = length(singleton_cost)
    keep_segment = falses(n)
    for i in 1:n
        if singleton_cost[i] < BIG_M
            keep_segment[i] = true
        end
    end

    feasible_segments = segments[keep_segment]
    feasible_singleton_cost = singleton_cost[keep_segment]
    feasible_singleton_boards = singleton_boards[keep_segment]
    infeasible = (1:length(singleton_cost))[singleton_cost .== BIG_M]

    original_indices = findall(keep_segment)

    cost_matrix_vertical_clearance, cost_matrix_max_height, cost_matrix_unused_area, best_boards =
        get_matching_cost_matrix(feasible_segments, board_depths=board_depths)

    return cost_matrix_vertical_clearance, cost_matrix_max_height, cost_matrix_unused_area, best_boards, feasible_singleton_cost, feasible_singleton_boards, infeasible, original_indices
end

function get_objective_function_lumber_matching(x_maxs, min_heights, moment, b, type, board_depths)

    function objective_function_lumber_matching(all_interior_knots)

        all_interior_knots = collect(all_interior_knots)

        knot_points = []
        n_interior_knots = length(all_interior_knots)/length(x_maxs)
        n_interior_knots = Int(n_interior_knots)

        for i in 1:lastindex(x_maxs)
            # Get the interior knots for this beam
            start_idx = (i-1)*n_interior_knots + 1
            end_idx = i*n_interior_knots
            interior_knots = all_interior_knots[start_idx:end_idx]
            
            # Add 0 at start and max x at end
            beam_knots = [0.0; interior_knots; x_maxs[i][end]]
            append!(knot_points, beam_knots)
        end

        # 1. Fit envelope and extract beam segments
        all_x_knots, all_y_knots, all_segment_lengths = get_lumber_segments(knot_points, x_maxs, min_heights, moment=moment, b=b, type=type)
        segments = get_individual_segments(all_x_knots, all_y_knots, all_segment_lengths)
        cost_matrix_vertical_clearance, cost_matrix_max_height, cost_matrix_unused_area, matches_boards, singleton_cost, singleton_boards, infeasible, index_map = get_costs(segments, board_depths=board_depths)
        println("GOT MATRICES")
        matches, singletons, cost = match_segments(cost_matrix_unused_area, singleton_cost)
        println("MATCHED SEGMENTS")
        println("COST: $cost")

        return cost  # This is your scalar objective value
    end

    return objective_function_lumber_matching

end
