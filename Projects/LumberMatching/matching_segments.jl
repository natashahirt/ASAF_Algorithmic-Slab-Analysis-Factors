function get_lumber_segments(knot_points, x_maxs, min_heights; moment=[], b=20, type=:C24)
    all_x_knots = []
    all_y_knots = []
    all_segment_lengths = []
    n_knots_per_beam = Int(length(knot_points) / length(x_maxs))
    for i in 1:lastindex(x_maxs)
        beam_knots = Float64.(knot_points[((i-1)*n_knots_per_beam + 1):(i*n_knots_per_beam)])
        if isempty(moment)
            x_knots, y_knots, segment_lengths = fit_segments(x_maxs[i], min_heights[i], beam_knots)
        else
            x_knots, y_knots, segment_lengths = fit_segments(x_maxs[i], min_heights[i], beam_knots, moment=moment[i], b=b, type=type)
        end
        push!(all_x_knots, x_knots)
        push!(all_y_knots, y_knots)
        push!(all_segment_lengths, segment_lengths)
    end
    return all_x_knots, all_y_knots, all_segment_lengths
end

function get_individual_segments(x_knots, y_knots, segment_lengths)
    segments = []
    for i in 1:lastindex(x_knots)
        for j in 1:lastindex(x_knots[i])-1
            segment_l = segment_lengths[i][j]
            min_y = min(y_knots[i][j], y_knots[i][j+1])
            max_y = max(y_knots[i][j], y_knots[i][j+1])
            if min_y < 1e-3 && max_y < 1e-3
                continue
            end
            push!(segments, (min_y, max_y, segment_l))
        end
    end
    return segments
end

function get_best_board(match::NTuple{2,Int}, board_depths::Vector{Float64}, cost_matrix_max_height::Matrix{Float64})
    best_board_idx = get_best_board(cost_matrix_max_height[match[1], match[2]], board_depths)
    return best_board_idx
end

function get_best_board(max_height::Float64, board_depths::Vector{Float64})
    best_board_idx = findfirst(board_depths .> max_height)
    if isnothing(best_board_idx)
        best_board = -1
    else
        best_board = board_depths[best_board_idx]
    end
    return best_board
end

"""
    get_matching_cost(s1::NTuple{3,Number}, s2::NTuple{3,Number}; board_depths::Vector{Real}=Real[])

Calculate the cost metrics for matching two beam segments.

# Arguments
- `s1::NTuple{3,Number}`: First segment parameters (start_height, end_height, length)
- `s2::NTuple{3,Number}`: Second segment parameters (start_height, end_height, length)
- `board_depths::Vector{Real}`: Available board depths (optional)

# Returns
- `Tuple{Float64,Float64,Float64}`: (vertical_clearance, max_height, unused_area)
"""
function polygon_area(x, y)
    area = 0.5 * abs(sum(x[i] * y[i+1] - x[i+1] * y[i] for i in 1:length(x) - 1))
    return area
end

function get_matching_cost(s1::Union{NTuple{3,Number}, NTuple{3,Any}}, s2::Union{NTuple{3,Number}, NTuple{3,Any}}; board_depths::Vector{Float64}=Float64[])
    h11, h12, l1 = Float64.(s1)
    h21, h22, l2 = Float64.(s2)

    m1 = (h12 - h11) / l1

    vertical_clearance = max(h11 + (m1 * l2) - abs(h21 - h22), h11)
    max_height = maximum([h12, h22 + vertical_clearance])

    unused_area_start = h11 + (m1 * l2)

    polygon_1_x = [l2, l2, l1, l1, l2]
    polygon_1_y = [unused_area_start, max_height, max_height, h12, unused_area_start]
    polygon_2_x = [0, 0, l2, l2, 0]
    polygon_2_y = [h11, vertical_clearance, max_height-h21, unused_area_start, h11]

    polygon_1_area = polygon_area(polygon_1_x, polygon_1_y)
    polygon_2_area = polygon_area(polygon_2_x, polygon_2_y)

    unused_area = polygon_1_area + polygon_2_area

    if !isempty(board_depths)
        best_board = get_best_board(max_height, board_depths)
        if best_board == -1
            unused_area = BIG_M
        else
            unused_rect_area = (best_board - max_height) * l1
            unused_area += unused_rect_area
        end
    else
        best_board = max_height
    end

    return vertical_clearance, max_height, unused_area, best_board
end

"""
    get_singleton_cost(segment::NTuple{3,Number})

Calculate the cost for a single beam segment.

# Arguments
- `segment::NTuple{3,Number}`: Segment parameters (start_height, end_height, length)

# Returns
- `Float64`: Cost of the singleton segment
"""
function get_singleton_cost(segment::Union{NTuple{3,Number}, NTuple{3,Any}}; board_depths::Vector{Float64}=Float64[])
    h1, h2, l = Float64.(segment)
    unused_area = 0.5 * (h2 - h1) * l
    if !isempty(board_depths)
        best_board = get_best_board(h2, board_depths)
        if best_board == -1
            unused_area = BIG_M
        else
            unused_area += (best_board - h2) * l
        end
    else
        best_board = h2
    end
    return unused_area, best_board
end

"""
    get_matching_cost_matrix(segments::Vector{NTuple{3,Number}}; board_depths::Vector{Real}=Real[])

Calculate cost matrices for all possible segment pairs.

# Arguments
- `segments::Vector{NTuple{3,Number}}`: Vector of segment parameters
- `board_depths::Vector{Real}`: Available board depths (optional)

# Returns
- `Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64}}`: (vertical_clearance_matrix, max_height_matrix, unused_area_matrix)
"""
function get_matching_cost_matrix(segments::Union{Vector{Any}, Vector{NTuple{3,Number}}}; board_depths::Vector{Float64}=Float64[])
    n = length(segments)
    cost_matrix_vertical_clearance = zeros(Float64, n, n)
    cost_matrix_max_height = zeros(Float64, n, n) 
    cost_matrix_unused_area = zeros(Float64, n, n)
    best_boards = zeros(Float64, n, n)

    if !isempty(board_depths)
        board_depths = sort(board_depths)
    end

    for i in 1:n
        for j in i+1:n # Only compute upper triangle
            cost_matrix_vertical_clearance[i,j], cost_matrix_max_height[i,j], cost_matrix_unused_area[i,j], best_boards[i,j] = get_matching_cost(segments[i], segments[j], board_depths=board_depths)
            # Mirror values to lower triangle
            cost_matrix_vertical_clearance[j,i] = cost_matrix_vertical_clearance[i,j]
            cost_matrix_max_height[j,i] = cost_matrix_max_height[i,j]
            cost_matrix_unused_area[j,i] = cost_matrix_unused_area[i,j]
            best_boards[j,i] = best_boards[i,j]
        end
    end

    return cost_matrix_vertical_clearance, cost_matrix_max_height, cost_matrix_unused_area, best_boards
end

"""
    get_singleton_cost_vector(segments::Vector{NTuple{3,Number}})

Calculate cost vector for all singleton segments.

# Arguments
- `segments::Vector{NTuple{3,Number}}`: Vector of segment parameters

# Returns
- `Vector{Float64}`: Vector of singleton costs
"""
function get_singleton_cost_vector(segments::Union{Vector{Any}, Vector{NTuple{3,Number}}}; board_depths::Vector{Float64}=Float64[])
    if !isempty(board_depths)
        board_depths = sort(board_depths)
    end
    costs_and_boards = get_singleton_cost.(segments, board_depths=board_depths)
    costs = [x[1] for x in costs_and_boards]
    boards = [x[2] for x in costs_and_boards]
    return costs, boards
end

