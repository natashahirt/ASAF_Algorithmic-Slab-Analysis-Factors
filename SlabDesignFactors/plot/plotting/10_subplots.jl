"""
    plot_10_subplots(df_all, subplot=:beam_sizer)

Creates a multi-plot figure to visualize various slab design factors.
"""
function plot_10_subplots(df_all; subplot::Symbol=:beam_sizer)

    @assert subplot in [:beam_sizer, :slab_sizer, :collinearity, :max_depth, :slab_type, :slab_min] "Invalid subplot type. Must be one of :beam_sizer, :slab_sizer, :collinearity, :max_depth, :slab_type"

    fig = Figure(size=(190*4, 190*4))

    # Common axis settings
    axis_kwargs = Dict(
        :aspect => DataAspect(),
        :limits => (0,150,0,150), #(0, maximum(df_all.steel_ec) * 1.25, 0, maximum(df_all.slab_ec) * 1.25),
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

    # Filter out rows with NaN total_ec
    df_all = filter(row -> !isnan(row.total_ec), df_all)

    # Create grids
    grid_main = fig[1,1] = GridLayout()

    function create_scatter(ax, df_usual, df_optimal; title="", ylabel=nothing, xlabel=nothing)
        # Plot scatter for usual and optimal data
        scatter!(ax, df_optimal.steel_ec, df_optimal.slab_ec, 
                 marker=df_optimal.symbol, rotation=df_optimal.rotation, 
                 color=(色[:irispurple]), transparency=sk[:transparency], 
                 markersize=sk[:markersize_large], alpha=sk[:alpha],
                 inspector_label=(self, i, p) -> df_optimal.category[i] * ": " * df_optimal.name[i])

        hull = andrew_hull(df_optimal.steel_ec, df_optimal.slab_ec)[1]
        points = [Point2f(df_optimal.steel_ec[i], df_optimal.slab_ec[i]) for i in hull]
        poly!(ax, points, color=(色[:irispurple], 0.05), transparency=true, strokewidth=0.5, strokecolor=色[:irispurple], linestyle=:dash)
    
        scatter!(ax, df_usual.steel_ec, df_usual.slab_ec, 
                 marker=df_usual.symbol, rotation=df_usual.rotation, 
                 color=(色[:skyblue]), transparency=sk[:transparency], 
                 markersize=sk[:markersize_large], alpha=sk[:alpha],
                 inspector_label=(self, i, p) -> df_usual.category[i] * ": " * df_usual.name[i])

        hull = andrew_hull(df_usual.steel_ec, df_usual.slab_ec)[1]
        points = [Point2f(df_usual.steel_ec[i], df_usual.slab_ec[i]) for i in hull]
        poly!(ax, points, color=(色[:skyblue], 0.05), transparency=true, strokewidth=0.5, strokecolor=色[:skyblue], linestyle=:dash)

        # Set axis properties
        ax.title = title
        if ylabel !== nothing
            ax.ylabel = ylabel
        end
        if xlabel !== nothing
            ax.xlabel = xlabel
        end
    end

    scatter_ax = Axis(grid_main[2, 1], tellheight=true, tellwidth=true; axis_kwargs...)

    if subplot != :slab_type
        if subplot == :beam_sizer
            df_usual = filter(row -> row.beam_sizer == "discrete", df_all)
            df_optimal = filter(row -> row.beam_sizer == "continuous", df_all)
            title = "Beam sizing"
            usual_label = "Catalog (W)"
            optimal_label = "Continuous"
        elseif subplot == :slab_sizer
            df_usual = filter(row -> row.slab_sizer == "uniform", df_all)
            df_optimal = filter(row -> row.slab_sizer == "cellular", df_all)
            title = "Slab sizing"
            usual_label = "Uniform"
            optimal_label = "Cellular"
        elseif subplot == :collinearity
            df_usual = filter(row -> row.collinear == true, df_all)
            df_optimal = filter(row -> row.collinear == false, df_all)
            title = "Beam collinearity"
            usual_label = "Collinear"
            optimal_label = "Noncollinear"
        elseif subplot == :max_depth
            df_usual = filter(row -> row.max_depth == 25, df_all)
            df_optimal = filter(row -> row.max_depth == 40, df_all)
            title = "Assembly depth"
            usual_label = "25\""
            optimal_label = "40\""
        elseif subplot == :slab_min
            df_usual = filter(row -> row.slab_min == true, df_all)
            df_optimal = filter(row -> row.slab_min == false, df_all)
            title = "Minimum slab depth"
            usual_label = "0.125m"
            optimal_label = "0.001m"
        end

        # Create axes and scatter plots
        create_scatter(scatter_ax, df_usual, df_optimal, xlabel="EC steel [kgCO2e/m²]", ylabel="EC RC-slab [kgCO2e/m²]")
        axislegend(scatter_ax, [sk[:elem_usual], sk[:elem_optimal]], [usual_label, optimal_label], position=:rt, orientation=:vertical, labelhalign=:left, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=sk[:fontsize], patchsize=(2, 10), padding=(0, 0, 0, 0))

        axis_x = Axis(grid_main[1, 1], limits=(0,150,nothing,nothing), height=25, title=title, titlesize=axis_kwargs[:titlesize])
        axis_y = Axis(grid_main[2, 2], limits=(nothing,nothing,0,150), width=25)

        values_x = [df_usual.steel_ec; df_optimal.steel_ec]
        categories_x = [[1 for _ in df_usual.steel_ec]; [2 for _ in df_optimal.steel_ec]]
        values_y = [df_usual.slab_ec; df_optimal.slab_ec]
        categories_y = [[1 for _ in df_usual.slab_ec]; [2 for _ in df_optimal.slab_ec]]
        
        color_x = map(d->d==1 ? 色[:skyblue] : 色[:irispurple], categories_x)
        color_y = map(d->d==1 ? 色[:skyblue] : 色[:irispurple], categories_y)

        boxplot!(axis_x, categories_x, values_x, orientation=:horizontal, color=color_x, show_notch=true, medianlinewidth=1, whiskerlinewidth=1, strokewidth=1, marker=:diamond, markersize=sk[:markersize_large], gap=0.5)
        boxplot!(axis_y, categories_y, values_y, color=color_y, show_notch=true, medianlinewidth=1, whiskerlinewidth=1, strokewidth=1, marker=:diamond, markersize=sk[:markersize_large], gap=0.5)

        hidedecorations!(axis_x)
        hidedecorations!(axis_y)
        hidespines!(axis_x)
        hidespines!(axis_y)
        axis_x.xgridvisible = true
        axis_y.ygridvisible = true

        steel_label_height = 2

    else        
        isotropic_data = filter(row -> row.slab_type == "isotropic", df_all)
        scatter!(scatter_ax, isotropic_data.steel_ec, isotropic_data.slab_ec, marker=:star8, color=(色[:skyblue], 0.2), transparency=true, markersize=sk[:markersize_large])
        orth_biaxial_data = filter(row -> row.slab_type == "orth_biaxial", df_all)
        scatter!(scatter_ax, orth_biaxial_data.steel_ec, orth_biaxial_data.slab_ec, marker=:cross, color=(色[:irispurple], 0.2), transparency=true, markersize=sk[:markersize_large])
        uniaxial_data = filter(row -> row.slab_type == "uniaxial", df_all)
        scatter!(scatter_ax, uniaxial_data.steel_ec, uniaxial_data.slab_ec, marker=:hline, color=(色[:magenta], 0.2), transparency=true, markersize=sk[:markersize_large], rotation=uniaxial_data.rotation)
        scatter_ax.xlabel = "EC steel [kgCO2e/m²]"
        scatter_ax.ylabel = "EC RC-slab [kgCO2e/m²]"

        hull = andrew_hull(orth_biaxial_data.steel_ec, orth_biaxial_data.slab_ec)[1]
        points = [Point2f(orth_biaxial_data.steel_ec[i], orth_biaxial_data.slab_ec[i]) for i in hull]
        poly!(scatter_ax, points, color=(色[:irispurple], 0.05), transparency=true, strokewidth=0.5, strokecolor=色[:irispurple], linestyle=:dash)

        hull = andrew_hull(uniaxial_data.steel_ec, uniaxial_data.slab_ec)[1]
        points = [Point2f(uniaxial_data.steel_ec[i], uniaxial_data.slab_ec[i]) for i in hull]
        poly!(scatter_ax, points, color=(色[:magenta], 0.02), transparency=true, strokewidth=0.5, strokecolor=色[:magenta], linestyle=:dash)

        hull = andrew_hull(isotropic_data.steel_ec, isotropic_data.slab_ec)[1]
        points = [Point2f(isotropic_data.steel_ec[i], isotropic_data.slab_ec[i]) for i in hull]
        poly!(scatter_ax, points, color=(色[:skyblue], 0.05), transparency=true, strokewidth=0.5, strokecolor=色[:skyblue], linestyle=:dash)

        axislegend(scatter_ax, [sk[:elem_isotropic], sk[:elem_orthogonal], sk[:elem_uniaxial]], ["Isotropic", "Biaxial Orthogonal", "Uniaxial"], position=:rt, orientation=:vertical, labelhalign=:left, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=sk[:fontsize], patchsize=(2, 10), padding=(0, 0, 0, 0))

        title = "Slab types"

        axis_x = Axis(grid_main[1, 1], limits=(0,150,nothing,nothing), height=25, title=title, titlesize=axis_kwargs[:titlesize])
        axis_y = Axis(grid_main[2, 2], limits=(nothing,nothing,0,150), width=25)

        values_x = [isotropic_data.steel_ec; orth_biaxial_data.steel_ec; uniaxial_data.steel_ec]
        categories_x = [[1 for _ in isotropic_data.steel_ec]; [2 for _ in orth_biaxial_data.steel_ec]; [3 for _ in uniaxial_data.steel_ec]]
        values_y = [isotropic_data.slab_ec; orth_biaxial_data.slab_ec; uniaxial_data.slab_ec]
        categories_y = [[1 for _ in isotropic_data.slab_ec]; [2 for _ in orth_biaxial_data.slab_ec]; [3 for _ in uniaxial_data.slab_ec]]
        
        color_x = map(d->d==1 ? 色[:skyblue] : d==2 ? 色[:irispurple] : 色[:magenta], categories_x)
        color_y = map(d->d==1 ? 色[:skyblue] : d==2 ? 色[:irispurple] : 色[:magenta], categories_y)

        boxplot!(axis_x, categories_x, values_x, orientation=:horizontal, color=color_x, show_notch=true, medianlinewidth=1, whiskerlinewidth=1, strokewidth=1, marker=:diamond, markersize=sk[:markersize_large], gap=0.5)
        boxplot!(axis_y, categories_y, values_y, color=color_y, show_notch=true, medianlinewidth=1, whiskerlinewidth=1, strokewidth=1, marker=:diamond, markersize=sk[:markersize_large], gap=0.5)

        hidedecorations!(axis_x)
        hidedecorations!(axis_y)
        hidespines!(axis_x)
        hidespines!(axis_y)
        axis_x.xgridvisible = true
        axis_y.ygridvisible = true

        steel_label_height = 1

    end

    # Find the business-as-usual (BAU) scenario
    bau_filter = row -> row.name == "r1c2" && row.slab_type == "uniaxial" && row.slab_sizer == "uniform" && row.beam_sizer == "discrete" && row.collinear == true && row.vector_1d_x == 1 && row.vector_1d_y == 0 && row.max_depth == 40
    business_as_usual = filter(bau_filter, df_all)

    bau_steel = business_as_usual.steel_ec
    bau_slab = business_as_usual.concrete_ec + business_as_usual.rebar_ec
    bau_total = bau_steel[1] + bau_slab[1]

    # Add BAU lines to all axes
    vlines!(scatter_ax, bau_steel, color=:black, linestyle=:dash, transparency=true, linewidth=1)
    hlines!(scatter_ax, bau_slab, color=:black, linestyle=:dash, transparency=true, linewidth=1)
    vlines!(axis_x, bau_steel, color=:black, linestyle=:dash, transparency=true, linewidth=1)
    hlines!(axis_y, bau_slab, color=:black, linestyle=:dash, transparency=true, linewidth=1)
    lines!(scatter_ax, [bau_total, 0], [0, bau_total], color=:black, linestyle=:dash, transparency=true, linewidth=1)
    text!(scatter_ax, 149, bau_slab[1] + 1, text="BAU slab", color = :black, align = (:right, :bottom), fontsize = sk[:fontsize_small])
    text!(scatter_ax, bau_steel[1] + 2, steel_label_height, text="BAU steel", color = :black, rotation = pi/2, align = (:left, :center), fontsize = sk[:fontsize_small])
    text!(scatter_ax, 2, bau_slab[1]+bau_steel[1] + 1, text="BAU total", color = :black, rotation=-pi/4, align = (:left, :center), fontsize = sk[:fontsize_small])    

    linkxaxes!(scatter_ax, axis_x)
    linkyaxes!(scatter_ax, axis_y)

    di = DataInspector(fig)

    GC.gc()

    return fig

end
