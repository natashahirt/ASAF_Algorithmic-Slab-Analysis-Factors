"""
    plot_5_barplot_maxdepth(df::DataFrame; category=nothing)

Plots a bar chart showing the distribution of beam sizes by section for specified maximum depths.
If a category is provided, the data is filtered to include only that category.

# Arguments
- `df::DataFrame`: The input data frame containing beam information.
- `category`: An optional category to filter the data by.
"""
function plot_5_beam_sizes(df_input::DataFrame; category=nothing)
    
    # Filter by category if provided
    if !isnothing(category)
        df_all = filter(row -> row.category == category, df_input)
    else
        df_all = df_input
    end

    # Plotting setup
    fontsize = 11
    smallfontsize = 8
    fig = Figure(size=(190*4, 190*2))
    grid = GridLayout(fig[2, 1])
    Label(grid[0, 1:2], text="b) Beam sizing by area", fontsize=fontsize, font=:bold)

    # Filter for discrete beam sizes
    df = filter(row -> row.beam_sizer == "discrete", df_all)
    df_names = DataFrame(name=String[], count=Int64[], percent=Float64[], depth=Float64[], mass=Float64[], max_depth=Int64[])

    # Section plot

    # Process for each max depth
    for max_depth in [25, 40]
        W_names = String[]

        # Collect names for the current max depth
        for row in eachrow(df)
            if row.max_depth == max_depth
                append!(W_names, parse_ids(row.ids))
            end
        end

        # Count occurrences and calculate mass
        dict_name_count = countmap(W_names)
        dict_name_mass = Dict(key => 0 for key in keys(dict_name_count))

        for W_name in keys(dict_name_count)
            split_name = split(W_name, r"(?<=\d)(?=\D)|(?<=\D)(?=\d)")
            W_depth = parse(Int, split_name[2])
            W_mass = parse(Int, split_name[4])
            W_count = dict_name_count[W_name]
            W_percent = (W_count / length(W_names) * 100) / 2
            push!(df_names, (W_name, W_count, W_percent, W_depth, W_mass, max_depth))
        end

        # Sort and filter for top 25% by percent
        sort!(df_names, [:percent], rev=true)
        top_50_percent_index = Int(ceil(0.5 * nrow(df_names)))
        df_names = df_names[1:70, :]

        sort!(df_names, [:mass, :depth])
    end

    # Assign categories to names
    seen = String[]
    categories = Int64[]
    category_dict = Dict()

    for name in df_names.name
        if name in seen
            push!(categories, findfirst(x -> x == name, seen))
        else
            new_category = isempty(categories) ? 1 : maximum(categories) + 1
            push!(categories, new_category)
            push!(seen, name)
            category_dict[new_category] = name
        end
    end

    df_names.category .= categories

    # Create axes for the plots
    ax = Axis(fig[1, 1], 
        title="a) Beam sizing by section", 
        xticks=(unique(df_names.category), [category_dict[category] for category in unique(df_names.category)]), 
        ylabel="% of total sections", 
        xticklabelrotation=-pi/4, 
        topspinevisible=false, 
        rightspinevisible=false, 
        yticklabelsize=fontsize, 
        xticklabelsize=smallfontsize, 
        xlabelsize=fontsize, 
        ylabelsize=fontsize, 
        titlesize=fontsize
    )

    # Plot grouped bar charts for each max depth
    max_depths = [25, 40]
    colors = [色[:skyblue], 色[:magenta]]

    for (i, max_depth) in enumerate(max_depths)
        color = colors[i]
        df_max_depth = filter(row -> row.max_depth == max_depth, df_names)
        barplot!(ax, df_max_depth.category .+ (i-1)*0.4, df_max_depth.percent, color=color, width=0.5, inspector_label=(self, j, p) -> df_max_depth.name[j], transparency=true)
    end

    # Legend setup
    elem_25 = MarkerElement(color=(色[:skyblue]), marker=:rect)
    elem_40 = MarkerElement(color=(色[:magenta]), marker=:rect)
    axislegend(ax, [elem_25, elem_40], ["25\"", "40\""], position=:rt, orientation=:vertical, labelhalign=:right, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=fontsize)

   # Area plot ================================

    # Initialize DataFrame for names, areas, and beam sizer types
    df_names = DataFrame(name=String[], area=Float64[], beam_sizer=String[])

    # Process each beam sizer type
    for beam_sizer in ["discrete", "continuous"]
        W_ids = [parse_ids(row.ids) for row in eachrow(df_all) if row.beam_sizer == beam_sizer]
        W_ids = vcat(W_ids...)  # Flatten the list of lists

        for W_id in W_ids
            if occursin("W", W_id)
                split_name = split(W_id, r"(?<=\d)(?=\D)|(?<=\D)(?=\d)") # add slashes before the D's
                W_area = W_imperial(W_id).A
                beam_sizer = "discrete"
            else
                W_area = parse(Float64, W_id)
                beam_sizer = "continuous"
            end
            push!(df_names, (W_id, W_area, beam_sizer))
        end
    end

    # Assign categories
    category_dict = Dict{String, Int}()
    df_names.category = [get!(category_dict, name, length(category_dict) + 1) for name in df_names.name]

    # Plotting

    # Plot discrete beam sizer
    ax1 = Axis(grid[1, 1], xticklabelrotation=pi/2, topspinevisible=false, rightspinevisible=false,
               xlabel="Area [in²]", ylabel="% of sections", limits=(0, 40, 0, nothing),
               yticks=(0:0.1:1, [string(label) for label in collect(0:10:100)]),
               yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize)
    df_filtered = filter(row -> row.beam_sizer == "discrete", df_names)
    hist!(ax1, df_filtered.area, color=色[:skyblue], normalization=:probability, bins=100)

    # Plot continuous beam sizer
    ax2 = Axis(grid[1, 2], xticklabelrotation=pi/2, topspinevisible=false, rightspinevisible=false,
               xlabel="Area [in²]", limits=(0, 40, 0, nothing),
               yticks=(0:0.1:1, [string(label) for label in collect(0:10:100)]),
               yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize)
    df_filtered = filter(row -> row.beam_sizer == "continuous", df_names)
    hist!(ax2, df_filtered.area, color=色[:magenta], normalization=:probability, bins=100)

    # Legend
    elem_discrete = MarkerElement(color=色[:skyblue], marker=:rect)
    elem_continuous = MarkerElement(color=色[:magenta], marker=:rect)
    axislegend(ax2, [elem_discrete, elem_continuous], ["Discrete", "Continuous"],
               position=:rt, orientation=:vertical, labelhalign=:right, framevisible=true,
               backgroundcolor=:white, framecolor=:white, labelsize=fontsize)

    # Data Inspector
    di = DataInspector(fig)

    # Link axes
    linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)

    GC.gc()

    return fig

end


function plot_5_beam_sizes_topology(df_input::DataFrame)
    
    # Filter by category if provided
    categories = [nothing, "nova", "topology", "grid"]
    titles = ["a) All layouts", "b) Nova", "c) Topology", "d) Grid"]
    tuples = [[1,1],[1,2],[2,1],[2,2]]
    axes = []

    # Plotting setup
    fontsize = 11
    smallfontsize = 8
    fig = Figure(size=(190*4, 190*2))
    grid = GridLayout(fig[1, 1])

    for (i,category) in enumerate(categories)

        if !isnothing(category)
            df_all = filter(row -> row.category == category, df_input)
        else
            df_all = df_input
        end

        # Filter for discrete beam sizes
        df = filter(row -> row.beam_sizer == "discrete", df_all)
        df_names = DataFrame(name=String[], count=Int64[], percent=Float64[], depth=Float64[], mass=Float64[], max_depth=Int64[])

        # Section plot

        # Process for each max depth
        for max_depth in [25, 40]
            W_names = String[]

            # Collect names for the current max depth
            for row in eachrow(df)
                if row.max_depth == max_depth
                    append!(W_names, parse_ids(row.ids))
                end
            end

            println(category, "...", length(W_names))

            # Count occurrences and calculate mass
            dict_name_count = countmap(W_names)

            for W_name in keys(dict_name_count)
                split_name = split(W_name, r"(?<=\d)(?=\D)|(?<=\D)(?=\d)")
                W_depth = parse(Int, split_name[2])
                W_mass = parse(Int, split_name[4])
                W_count = dict_name_count[W_name]
                W_percent = (W_count / length(W_names) * 100) / 2
                push!(df_names, (W_name, W_count, W_percent, W_depth, W_mass, max_depth))
            end

            # Sort and filter for top 25% by percent
            sort!(df_names, [:percent], rev=true)
            #top_50_percent_index = Int(ceil(0.5 * nrow(df_names)))
            df_names = df_names[1:35, :]

            sort!(df_names, [:mass, :depth])
        end

        # Assign categories to names
        seen = String[]
        categories = Int64[]
        category_dict = Dict()

        for name in df_names.name
            if name in seen
                push!(categories, findfirst(x -> x == name, seen))
            else
                new_category = isempty(categories) ? 1 : maximum(categories) + 1
                push!(categories, new_category)
                push!(seen, name)
                category_dict[new_category] = name
            end
        end

        df_names.category .= categories

        # Create axes for the plots
        ax = Axis(grid[tuples[i][1], tuples[i][2]], 
            title=titles[i],
            xticks=(unique(df_names.category), [category_dict[category] for category in unique(df_names.category)]), 
            ylabel="% of total sections", 
            xticklabelrotation=-pi/4, 
            topspinevisible=false, 
            rightspinevisible=false, 
            yticklabelsize=fontsize, 
            xticklabelsize=smallfontsize, 
            xlabelsize=fontsize, 
            ylabelsize=fontsize, 
            titlesize=fontsize
        )

        push!(axes, ax)

        # Plot grouped bar charts for each max depth
        max_depths = [25, 40]
        colors = [色[:skyblue], 色[:magenta]]

        for (i, max_depth) in enumerate(max_depths)
            color = colors[i]
            df_max_depth = filter(row -> row.max_depth == max_depth, df_names)
            barplot!(ax, df_max_depth.category .+ (i-1)*0.2, df_max_depth.percent, color=color, width=0.3, inspector_label=(self, j, p) -> df_max_depth.name[j], transparency=true)
        end

        # Legend setup
        if length(axes) == 2
            elem_25 = MarkerElement(color=(色[:skyblue]), marker=:rect)
            elem_40 = MarkerElement(color=(色[:magenta]), marker=:rect)
            axislegend(ax, [elem_25, elem_40], ["25\"", "40\""], position=:rt, orientation=:vertical, labelhalign=:right, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=fontsize)
        end

    end

    linkyaxes!(axes[1], axes[2])
    linkyaxes!(axes[3], axes[4])

    GC.gc()

    return fig

end
