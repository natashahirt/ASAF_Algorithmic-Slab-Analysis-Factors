using Random

for i in 1:20
    Random.seed!(i)  # Set random seed for reproducibility
    segment_idx_1 = rand(1:length(segments))
    segment_idx_2 = rand(1:length(segments))

    segment_1 = segments[segment_idx_1]
    segment_2 = segments[segment_idx_2]

    if segment_1[3] < segment_2[3] # longer segment on the bottom
        segment_1, segment_2 = segment_2, segment_1 
    end

    # Create a new figure for segment comparison
    fig2 = Figure();
    ax2 = Axis(fig2[1,1],
        xlabel = "Position along beam (in)",
        ylabel = "Height (in)", 
        title = "Comparison of Segments $segment_idx_1 and $segment_idx_2")

    # Plot segment 1
    x_1 = [0, 0, segment_1[3], segment_1[3], 0]  # x coordinates from 0 to length
    y_1 = [0, segment_1[1], segment_1[2], 0, 0]  # y coordinates from min to max height
    lines!(ax2, x_1, y_1, label="Segment $segment_idx_1", color=:blue)

    # Plot segment 2
    vertical_clearance, max_height, unused_area, best_board = get_matching_cost(segment_1, segment_2)
    println(max_height, " ", unused_area)
    x_2 = [0, 0, segment_2[3], segment_2[3], 0]
    y_2 = [vertical_clearance, max_height, max_height, max_height - segment_2[1], vertical_clearance]
    lines!(ax2, x_2, y_2, label="Segment $segment_idx_2", color=:red)

    display(fig2)
end


fig3 = Figure(size=(1500, 500));
ax3 = Axis(fig3[1,1],
    xlabel = "Segment 1",
    ylabel = "Segment 2",
    title = "Vertical Clearance",)

    ax4 = Axis(fig3[1,2],
    xlabel = "Segment 1",
    ylabel = "Segment 2",
    title = "Max Height",)
ax5 = Axis(fig3[1,3],
    xlabel = "Segment 1",
    ylabel = "Segment 2",
    title = "Unused Area",)

heatmap!(ax3, cost_matrix_vertical_clearance, label="Vertical Clearance", colormap=:magma)
heatmap!(ax4, cost_matrix_max_height, label="Max Height", colormap=:magma)
heatmap!(ax5, cost_matrix_unused_area, label="Unused Area", colormap=:magma)

display(fig3)


fig = Figure(size=(1000, 1000));
ax = Axis(fig[1,1], title="Segments", xlabel="Length (in)", ylabel="Depth (in)");
for i in 1:lastindex(all_x_knots)
    x = all_x_knots[i]
    y = all_y_knots[i]
    lines!(ax, x_maxs[i], min_heights[i], label="Beam $i", linestyle=:dash)
    lines!(ax, x, y, label="Beam $i", color=:grey)
end
fig
