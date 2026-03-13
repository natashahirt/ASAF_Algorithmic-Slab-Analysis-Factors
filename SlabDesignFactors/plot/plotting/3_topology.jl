"""
    plot_3_topology(df)

Creates a topology plot for slab design factors, visualizing the best and worst scenarios.
"""
function plot_3_topology(df; category=nothing)

    maximum_ec = maximum(filter(x -> !isnan(x), df.total_ec)) * 1.05

    if !isnothing(category)
        df = filter(row -> row.category == category, df)
    end

    # Prepare data
    df.shorthand_best .= ""
    df.shorthand_worst .= ""
    
    slab_names = unique(df.name)
    df_best = similar(df, 0)
    df_worst = similar(df, 0)

    for name in slab_names
        df_name = filter(row -> row.name == name, df)
        best_row = eachrow(df_name)[findmin(df_name.total_ec)[2]]
        worst_row = eachrow(df_name)[findmax(df_name.total_ec)[2]]

        # Best shorthand
        best_row.shorthand_best *= best_row.slab_sizer == "uniform" ? "u" : "c"
        best_row.shorthand_best *= best_row.beam_sizer == "discrete" ? "W" : "c"
        best_row.shorthand_best *= best_row.collinear ? "c" : "n"
        best_row.shorthand_best *= string(Int(best_row.max_depth))

        # Worst shorthand
        best_row.shorthand_worst *= worst_row.slab_sizer == "uniform" ? "u" : "c"
        best_row.shorthand_worst *= worst_row.beam_sizer == "discrete" ? "W" : "c"
        best_row.shorthand_worst *= worst_row.collinear ? "c" : "n"
        best_row.shorthand_worst *= string(Int(worst_row.max_depth))

        push!(df_best, best_row)
        push!(df_worst, worst_row)
    end

    df_best.worst_total_ec .= df_worst.total_ec
    df_worst.best_total_ec .= df_best.total_ec

    # Bar plot
    sort!(df_best, order(:total_ec))
    sort!(df_worst, order(:best_total_ec))

    if !isnothing(category)
        df_best.rowcol = replace.(df_best.rowcol, category => "t")
        df_worst.rowcol = replace.(df_worst.rowcol, category => "t")
    end

    # Create figure
    fig = Figure(size=(190*4, 190*2.5))
    grid = fig[1,1] = GridLayout()
    ax1 = Axis(grid[1,1], ylabel="EC [kgCO2e/m²]", xticklabelrotation=pi/2, 
               limits=(0, length(df_best.name) + 1, 0, maximum_ec), 
               xticks=(1:lastindex(df_best.name), df_best.rowcol), 
               yticklabelsize=11, xticklabelsize=11, xlabelsize=11, ylabelsize=11, titlesize=11)

    grid_topologies = grid[2,1] = GridLayout()
    rowsize!(grid, 2, Fixed(190*0.5))

    positions = 1:lastindex(df_best.area)
    stack = repeat(positions, inner=3)
    categories = repeat(positions, outer=3)
    colours = repeat([1, 2, 3], inner=lastindex(df_best.area))

    hatchpattern = Makie.LinePattern(direction=[1, 1]; width=2, tilesize=(5, 5), linecolor=:lightgrey, background_color=:white)
    barplot!(ax1, positions, df_best.worst_total_ec, color=hatchpattern, direction=:y, 
             bar_labels=df_best.shorthand_worst, label_font=:italic, label_size=7, label_rotation=pi/2,
             flip_labels_at=(maximum_ec * 0.8))

    df_concrete = df_best.concrete_ec
    df_steel = df_best.steel_ec
    df_rebar = df_best.rebar_ec

    data_barplot = vcat(df_concrete, df_steel, df_rebar)

    barplot!(ax1, positions, df_best.total_ec, color=:white, direction=:y, 
             bar_labels=df_best.shorthand_best, label_font=:italic, label_size=8, label_rotation=pi/2)
    barplot!(ax1, categories, data_barplot, stack=stack, color=colours, direction=:y, 
             colormap=[(色[:charcoalgrey]), (色[:skyblue]), (色[:magenta])], strokewidth=0.5)

    # let's try draw the best and worst topologies
    best_topologies = df_best[1:3,:]
    worst_topologies = df_best[length(df_best.name) - 2:end,:]

    println(worst_topologies)

    topology_names = vcat(best_topologies.name, worst_topologies.name)
    topology_rowcol = vcat(best_topologies.rowcol, worst_topologies.rowcol)

    for i in 1:6
        ax_topologies = Axis(grid_topologies[1,i], aspect=DataAspect())

        path = "Geometries/$category/$(topology_names[i]).json"
        geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))
        geometry = generate_from_json(geometry_dict, plot=false, drawn=true);

        plot_elements!(ax_topologies, geometry.elements, linewidth=1)
        
        if !isnothing(category)
            xlabel = replace(topology_rowcol[i], category => "t")
        else
            xlabel = topology_names[i]
        end

        hidedecorations!(ax_topologies, label=false)
        hidespines!(ax_topologies)

        ax_topologies.xlabel = xlabel
        ax_topologies.xlabelsize = 11
    end
    
    GC.gc()

    return fig

end
