function plot_mass_against_unique_sections(df_all)

    fig = Figure(size=(190*4, 190*1.5))

    # Common axis settings
    axis_kwargs = Dict(
        :yticklabelsize => 11,
        :xticklabelsize => 11,
        :xlabelsize => 11,
        :ylabelsize => 11,
        :titlesize => 11
    )

    sk = Dict(
        :transparency => true,
        :markersize => 1.5,
        :markersize_large => 6,
        :fontsize => 11,
        :fontsize_small => 8,
        :alpha => 0.4,
        :elem_usual => MarkerElement(color = 色[:skyblue], marker = :circle),
        :elem_optimal => MarkerElement(color = 色[:irispurple], marker = :circle),
        :elem_isotropic => MarkerElement(color = 色[:skyblue], marker = :star8),
        :elem_orthogonal => MarkerElement(color = 色[:irispurple], marker = :cross),
        :elem_uniaxial => MarkerElement(color = 色[:magenta], marker = :hline),
    )

    name_dict = Dict(
        "office" => "Office",
        "school" => "School",
        "warehouse" => "Warehouse",
    )

    x_lims_dict = Dict(
        "office" => (0, 21),
        "school" => (0, 13),
        "warehouse" => (0, 22),
    )

    y_lims_dict = Dict(
        "office" => (0, 150),
        "school" => (0, 300),
        "warehouse" => (0, 140),
    )

    # Filter out rows with NaN total_ec
    df_all = filter(row -> !isnan(row.total_ec), df_all)

    axes = Axis[]
    color_list = [色[:lilac], 色[:skyblue], 色[:magenta]]

    # Iterate over each category and create a subplot
    for (i, (key, display_name)) in enumerate(name_dict)
        
        df = filter(row -> row.name == key && row.unique_sections <= x_lims_dict[key][2], df_all)
        min_steel_norm = minimum(df.steel_norm)

        ax = Axis(fig[1, i]; axis_kwargs...)
        ylims!(ax, y_lims_dict[key]...)
        xlims!(ax, x_lims_dict[key]...)
        push!(axes, ax)  # Collect axes for linking

        # Calculate percentage increases
        pct_weight_increase = 100 .* (df.steel_norm ./ min_steel_norm .- 1)
        # Filter out zero values
        nonzero_mask = df.unique_sections .!= 0
        ln = lines!(ax, df.unique_sections[nonzero_mask], pct_weight_increase[nonzero_mask],
               color=color_list[i], transparency=sk[:transparency])
        pts = scatter!(ax, [NaN, NaN], marker=:circle, color=:black, markersize=5)

        axislegend(ax, [ln, pts], ["Novel method", "ETABS"], position=:rt, framevisible=true, framecolor=:white, labelsize=11)

        # Set custom ticks and labels
        max_val = x_lims_dict[key][2]
        ticks = [0, 5, 10, 15, 20]  # Start with standard increments
        # Only include values up to the max, then add the max itself
        ticks = [t for t in ticks if t < max_val]
        push!(ticks, max_val)
        ax.xticks = ticks
        ax.xtickformat = xs -> [t == max_val ? "∞" : string(t) for t in ticks]

        ax.xlabel = "# unique sections"
        if i == 1
            ax.ylabel = "Weight increase in steel grillage [%]"
        end
        ax.title = display_name
    end

    return fig
end

function plot_section_distribution(df_all; name="school", max_height="false")

    # Filter out rows with NaN total_ec
    df_all = filter(row -> !isnan(row.total_ec), df_all)

    #df = filter(row -> row.name == name && row.unique_sections == 30, df_all)
    df = filter(row -> row.name == name, df_all)
    
    csv_source = "Geometries/validation/$(name)_original.csv"
    df_csv = CSV.read(csv_source, DataFrame, header=false)
    # Add column headers if not present
    rename!(df_csv, [:depth, :weight])
    original_depths = df_csv.depth
    original_weights = df_csv.weight

    fig = Figure(size=(190*4, 190*1.5))

    # Common axis settings
    axis_kwargs = Dict(
        :yticklabelsize => 11,
        :xticklabelsize => 11,
        :xlabelsize => 11,
        :ylabelsize => 11,
        :titlesize => 11
    )

    # Create two subplots
    ax1 = Axis(fig[1, 1]; axis_kwargs...)
    ax2 = Axis(fig[1, 2]; axis_kwargs...)
    ylims!(ax1, (0, nothing))
    ylims!(ax2, (0, nothing))

    # Extract depths and weights from section names
    function get_section_data(section::String)
        m = match(r"W(\d+)X(\d+\.?\d*)", section)
        if m !== nothing
            return parse(Float64, m[1]), parse(Float64, m[2])
        end
        return 0.0, 0.0
    end

    # Get all sections and parse them
    all_sections = parse_sections(df.sections[1])
    depths = Float64[]
    weights = Float64[]
    for section in all_sections
        depth, weight = get_section_data(section)
        push!(depths, depth)
        push!(weights, weight)
    end

    # Define colors for each project type
    color_list = Dict(
        "school" => 色[:skyblue],
        "office" => 色[:lilac], 
        "warehouse" => 色[:magenta]
    )

    # Get color based on project name
    hist_color = color_list[name]

    # Compute the min and max across both datasets
    all_weights = vcat(weights, original_weights)
    bin_edges = range(minimum(all_weights), stop=maximum(all_weights), length=21)  # 20 bins

    # Use the same edges for both histograms, and assign to variables for legend
    h1 = hist!(ax1, weights, bins=bin_edges, color=hist_color, strokewidth=1, label="Novel method")
    h2 = hist!(ax1, original_weights, bins=bin_edges, color=(:black, 0.25), strokewidth=1, strokecolor=:black, label="ETABS")

    all_depths = vcat(depths, original_depths)
    bin_edges_depth = range(minimum(all_depths), stop=maximum(all_depths), length=21)  # 20 bins

    h3 = hist!(ax2, depths, bins=bin_edges_depth, color=hist_color, strokewidth=1, label="Novel method")
    h4 = hist!(ax2, original_depths, bins=bin_edges_depth, color=(:black, 0.25), strokewidth=1, strokecolor=:black, label="ETABS")

    # Labels and titles
    ax1.xlabel = "Linear Weight (lb/ft)"
    ax1.ylabel = "Number of Instances"
    ax1.title = "Distribution of Section Weights"

    ax2.xlabel = "Approximate Section Depth"
    ax2.ylabel = "Number of Instances"
    ax2.title = "Distribution of Section Depths"

    if max_height
        if df.max_depth[1] == 1000
            Label(fig[0, :], "MAX DEPTH = ∞ IN", fontsize=11, font=:bold)
        else
            Label(fig[0, :], "MAX DEPTH = $(df.max_depth[1]) IN", fontsize=11, font=:bold)
        end
    end

    # Add legend to the figure (for ax1, but you can also add for ax2 if you want)
    axislegend(ax2, [h1, h2], ["Novel method", "ETABS"], framecolor=:white, position=:lt, orientation=:vertical, labelsize=9, patchsize=(10.0f0, 10.0f0))

    return fig
end


function plot_mass_against_max_depth(df_all)

    fig = Figure(size=(190*4, 190*1.5))

    # Common axis settings
    axis_kwargs = Dict(
        :yticklabelsize => 11,
        :xticklabelsize => 11,
        :xlabelsize => 11,
        :ylabelsize => 11,
        :titlesize => 11
    )

    sk = Dict(
        :transparency => true,
        :markersize => 1.5,
        :markersize_large => 6,
        :fontsize => 11,
        :fontsize_small => 8,
        :alpha => 0.4,
        :elem_usual => MarkerElement(color = 色[:skyblue], marker = :circle),
        :elem_optimal => MarkerElement(color = 色[:irispurple], marker = :circle),
        :elem_isotropic => MarkerElement(color = 色[:skyblue], marker = :star8),
        :elem_orthogonal => MarkerElement(color = 色[:irispurple], marker = :cross),
        :elem_uniaxial => MarkerElement(color = 色[:magenta], marker = :hline),
    )

    name_dict = Dict(
        "office" => "Office",
        "school" => "School",
        "warehouse" => "Warehouse",
    )

    # Filter out rows with NaN total_ec
    df_all = filter(row -> !isnan(row.total_ec), df_all)

    axes = Axis[]
    color_list = [色[:lilac], 色[:skyblue], 色[:magenta]]

    # Iterate over each category and create a subplot
    for (i, (key, display_name)) in enumerate(name_dict)

        source_path = "Geometries/validation/$(display_name)_original.csv"
        df = CSV.read(source_path, DataFrame, header=false)
        rename!(df, [:depth, :weight])
        max_depth = maximum(df.depth) + 4.75

        println(display_name, " ", maximum(df.depth), " ", max_depth)
        
        df = filter(row -> row.name == key, df_all)
        
        ax = Axis(fig[1, i]; axis_kwargs...)
        xlims!(ax, (0, 55))
        ylims!(ax, (0, nothing))
        push!(axes, ax)  # Collect axes for linking

        # Calculate percentage increases
        scatterlines!(ax, df.max_depth, df.steel_norm,
               color=color_list[i], transparency=sk[:transparency])

        # Add vertical line at max depth
        vlines!(ax, max_depth, color=:black, linestyle=:dash, linewidth=1)

        # Set custom ticks and labels
        ticks = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]  # Start with standard increments
        ax.xticks = ticks
        ax.xtickformat = xs -> [t == 55 ? "∞" : string(t) for t in ticks]

        ax.xlabel = "Maximum assembly depth [in]"
        ax.ylabel = "Normalized mass of steel beams [kg/m²]"
        ax.title = display_name
    end

    linkaxes!(axes...)

    return fig
end