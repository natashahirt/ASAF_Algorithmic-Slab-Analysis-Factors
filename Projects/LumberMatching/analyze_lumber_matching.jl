include("../../SlabDesignFactors/scripts/_scripts.jl")

path = "Geometries/grid/x5y0.json"
path = "Geometries/topology/r9c4.json"

using NLopt
const GUROBI_ENV = Gurobi.Env()

name = basename(splitext(path)[1])    # Name for the plot
# Parse geometry from JSON
geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))
geometry, type_information = generate_from_json(geometry_dict, plot=true, drawn=false);

# Analyze the slab to get dimensions
slab_params = SlabAnalysisParams(
    geometry,
    slab_name=name,
    slab_type=:uniaxial,
    vector_1d=[1,0],
    slab_sizer=:uniform,
    spacing=.1,
    plot_analysis=true,
    fix_param=true,
    slab_units=:m,
);

# Sizing parameters
beam_sizing_params = SlabSizingParams(
    live_load=psf_to_ksi(50), # ksi
    superimposed_dead_load=psf_to_ksi(15), # ksi
    live_factor=1.6, # -
    dead_factor=1.2, # -
    beam_sizer=:discrete,
    max_depth=40, # in
    beam_units=:in, # in, etc.
    serviceability_lim=360,
    collinear=true,
    minimum_continuous=true,
    n_max_sections=0,
    slab_dead_load=0,
);

b = 20
type = :C24

M_maxs, V_maxs, x_maxs = get_internal_forces_slab(slab_params, beam_sizing_params); # only do this once
min_heights = get_min_heights(M_maxs, V_maxs, b=b, type=type); # only do this once

n_knots = 5
knot_points = []
for i in 1:length(x_maxs)
    beam_length = x_maxs[i][end]
    beam_knots = range(0, beam_length, length=n_knots)
    append!(knot_points, beam_knots)
end

all_x_knots, all_y_knots, all_segment_lengths = get_lumber_segments(knot_points, x_maxs, min_heights, moment=M_maxs, b=b, type=type);
segments = get_individual_segments(all_x_knots, all_y_knots, all_segment_lengths)
cost_matrix_vertical_clearance, cost_matrix_max_height, cost_matrix_unused_area, best_boards = get_matching_cost_matrix(segments)
single_cost = get_singleton_cost_vector(segments)

board_depths = sort(Float64[6, 8, 10, 12, 15])

cost_matrix_vertical_clearance, cost_matrix_max_height, cost_matrix_unused_area, matches_boards, singleton_cost, singleton_boards, infeasible, index_map = get_costs(segments, board_depths=board_depths)
matches, singletons, cost = match_segments(cost_matrix_unused_area, singleton_cost)

n_knots = 5  # or however many you want to vary

initial_knots = []
lb = []
ub = []
for i in 1:length(x_maxs)
    beam_length = x_maxs[i][end]
    beam_knots = range(0, beam_length, length=n_knots)
    append!(initial_knots, beam_knots[2:end-1])
    append!(lb, fill(0.0, n_knots-2))
    append!(ub, fill(x_maxs[i][end], n_knots-2))
end

length(initial_knots) == length(lb) == length(ub)

board_depths = sort(Float64[6, 8, 10, 12, 15])
cost_function = get_objective_function_lumber_matching(x_maxs, min_heights, M_maxs, b, type, board_depths)

model = Nonconvex.Model(cost_function)
Nonconvex.addvar!(model, lb, ub)

alg = NLoptAlg(:LN_COBYLA);
options = NLoptOptions(maxeval=1000, xtol_abs=1, ftol_abs=1)

r = Nonconvex.optimize(model, alg, initial_knots, options = options);

r.minimizer
r.minimum

knot_points = []

n_interior_knots = length(r.minimizer)/length(x_maxs)
n_interior_knots = Int(n_interior_knots)

for i in 1:lastindex(x_maxs)
    # Get the interior knots for this beam
    start_idx = (i-1)*n_interior_knots + 1
    end_idx = i*n_interior_knots
    interior_knots = r.minimizer[start_idx:end_idx]
    
    # Add 0 at start and max x at end
    beam_knots = [0.0; interior_knots; x_maxs[i][end]]
    append!(knot_points, beam_knots)
end

all_x_knots, all_y_knots, all_segment_lengths = get_lumber_segments(knot_points, x_maxs, min_heights, moment=M_maxs, b=b, type=type);
segments = get_individual_segments(all_x_knots, all_y_knots, all_segment_lengths)
cost_matrix_vertical_clearance, cost_matrix_max_height, cost_matrix_unused_area, best_boards = get_matching_cost_matrix(segments)
single_cost = get_singleton_cost_vector(segments)

board_depths = sort(Float64[6, 8, 10, 12, 15])

cost_matrix_vertical_clearance, cost_matrix_max_height, cost_matrix_unused_area, matches_boards, singleton_cost, singleton_boards, infeasible, index_map = get_costs(segments, board_depths=board_depths)
matches, singletons, cost = match_segments(cost_matrix_unused_area, singleton_cost)

begin
    grouped_matches = Dict{Float64, Vector{Tuple{Int,Int}}}()
    for match in matches
        board = matches_boards[match[1], match[2]]
        if !haskey(grouped_matches, board)
            grouped_matches[board] = Tuple{Int,Int}[]
        end
        push!(grouped_matches[board], match)
    end
    for single in singletons
        board = singleton_boards[single]
        if !haskey(grouped_matches, board)
            grouped_matches[board] = Tuple{Int,Int}[]
        end
        push!(grouped_matches[board], (single, -1))
    end
    
    if !isempty(infeasible)
        fig = Figure(size = (1000, (length(keys(grouped_matches)) + 1) * 200));
    else
        fig = Figure(size = (1000, length(keys(grouped_matches)) * 200));
    end;

    # Add title to figure
    Label(fig[0,:], "Lumber Matching Results Cost $(round(cost, digits=2)) in²", fontsize=12, tellwidth=false)

    row = 1
    row_axes = []

    for (depth, pairlist) in sort(collect(grouped_matches))
        ax = Axis(fig[row, 1], title = "Board Depth $(depth) in", xlabel = "Length (in)", ylabel = "Depth (in)", height = 100)
        current_x = 0.0

        for (i1, i2) in pairlist
            if i2 == -1
                s1 = segments[index_map[i1]]
                x1 = [current_x, current_x, current_x + s1[3], current_x + s1[3], current_x]
                y1 = [0, s1[1], s1[2], 0, 0]
                lines!(ax, x1, y1, color = :grey, label = "Seg $i1")
                current_x += s1[3]
            else
                s1 = segments[index_map[i1]]
                s2 = segments[index_map[i2]]

                # Sort segments so s1 is the longer one
                if s2[3] > s1[3]
                    s1, s2 = s2, s1
                    i1, i2 = i2, i1
                end

                clearance, max_height, unused_area, best_board = get_matching_cost(s1, s2)
                if max_height > depth
                    println(s1, " ", s2)
                end

                # Segment 1 (blue)
                x1 = [current_x, current_x, current_x + s1[3], current_x + s1[3], current_x]
                y1 = [0, s1[1], s1[2], 0, 0]
                lines!(ax, x1, y1, color = :blue, label = "Seg $i1")

                # Segment 2 (red), shifted up
                x2 = [current_x, current_x, current_x + s2[3], current_x + s2[3], current_x]
                y2 = [clearance, max_height, max_height, max_height - s2[1], clearance]
                lines!(ax, x2, y2, color = :red, label = "Seg $i2")

                current_x += max(s1[3], s2[3])
            end
        end
        hlines!(ax, depth, color=:grey, linestyle=:dash)
        xlims!(ax, (0, nothing))
        ylims!(ax, (0, nothing))
        push!(row_axes, ax)
        row += 1
    end

    # Add row for infeasible segments
    if !isempty(infeasible)
        ax = Axis(fig[row, 1], title = "Infeasible Segments", xlabel = "Length (in)", ylabel = "Depth (in)", height = 100)
        current_x = 0.0
        
        for i in infeasible
            s = segments[i]
            x = [current_x, current_x, current_x + s[3], current_x + s[3], current_x]
            y = [0, s[1], s[2], 0, 0]
            lines!(ax, x, y, color = :gray, label = "Seg $i")
            current_x += s[3]
        end
        
        xlims!(ax, (0, nothing))
        ylims!(ax, (0, nothing))
        push!(row_axes, ax)
    end

    linkyaxes!(row_axes...)
    linkxaxes!(row_axes...)

    display(fig)
end