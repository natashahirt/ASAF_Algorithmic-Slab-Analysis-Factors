"""
    plot_1_multiplot(df_all)

Creates a multi-plot figure to visualize various slab design factors.
"""
function plot_1_multiplot(df_all)

    fig = Figure(size=(190*4, 190*2))

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
        :markersize_large => 4,
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

    grid_subplots = grid_main[1, 1] = GridLayout()
    grid_slabtype = grid_main[1, 2] = GridLayout()

    colgap!(grid_subplots, 2)
    rowgap!(grid_subplots, 2)
    colgap!(grid_main, 2)
    rowgap!(grid_main, 2)

    function create_scatter(ax, df_usual, df_optimal; title, ylabel=nothing, xlabel=nothing)
        # Plot scatter for usual and optimal data
        scatter!(ax, df_optimal.steel_ec, df_optimal.slab_ec, 
                 marker=df_optimal.symbol, rotation=df_optimal.rotation, 
                 color=(色[:irispurple]), transparency=sk[:transparency], 
                 markersize=sk[:markersize], alpha=sk[:alpha],
                 inspector_label=(self, i, p) -> df_optimal.category[i] * ": " * df_optimal.name[i])

        hull = andrew_hull(df_optimal.steel_ec, df_optimal.slab_ec)[1]
        points = [Point2f(df_optimal.steel_ec[i], df_optimal.slab_ec[i]) for i in hull]
        poly!(ax, points, color=(色[:irispurple], 0.05), transparency=true, strokewidth=0.5, strokecolor=色[:irispurple], linestyle=:dash)
    
        scatter!(ax, df_usual.steel_ec, df_usual.slab_ec, 
                 marker=df_usual.symbol, rotation=df_usual.rotation, 
                 color=(色[:skyblue]), transparency=sk[:transparency], 
                 markersize=sk[:markersize], alpha=sk[:alpha],
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

    # Create axes and scatter plots
    ax1 = Axis(grid_subplots[1, 1]; axis_kwargs...)
    create_scatter(ax1, filter(row -> row.slab_sizer == "uniform", df_all), filter(row -> row.slab_sizer == "cellular", df_all), title=" a) Slab sizing", ylabel="EC RC-slab [kgCO2e/m²]")

    ax2 = Axis(grid_subplots[1, 2]; axis_kwargs...)
    create_scatter(ax2, filter(row -> row.beam_sizer == "discrete", df_all), filter(row -> row.beam_sizer == "continuous", df_all), title="b) Beam sizing")

    ax3 = Axis(grid_subplots[2, 1]; axis_kwargs...)
    create_scatter(ax3, filter(row -> row.collinear == true, df_all), filter(row -> row.collinear == false, df_all), title="c) Beam collinearity", xlabel="EC steel [kgCO2e/m²]", ylabel="EC RC-slab [kgCO2e/m²]")

    ax4 = Axis(grid_subplots[2, 2]; axis_kwargs...)
    create_scatter(ax4, filter(row -> row.slab_min == true, df_all), filter(row -> row.slab_min == false, df_all), title="d) Minimum slab depth", xlabel="EC steel [kgCO2e/m²]")

    ax5 = Axis(grid_slabtype[1, 1]; axis_kwargs...)
    isotropic_data = filter(row -> row.slab_type == "isotropic", df_all)
    scatter!(ax5, isotropic_data.steel_ec, isotropic_data.slab_ec, marker=:star8, color=(色[:skyblue], 0.2), transparency=true, markersize=sk[:markersize_large])
    orth_biaxial_data = filter(row -> row.slab_type == "orth_biaxial", df_all)
    scatter!(ax5, orth_biaxial_data.steel_ec, orth_biaxial_data.slab_ec, marker=:cross, color=(色[:irispurple], 0.2), transparency=true, markersize=sk[:markersize_large])
    uniaxial_data = filter(row -> row.slab_type == "uniaxial", df_all)
    scatter!(ax5, uniaxial_data.steel_ec, uniaxial_data.slab_ec, marker=:hline, color=(色[:magenta], 0.2), transparency=true, markersize=sk[:markersize_large], rotation=uniaxial_data.rotation)
    ax5.title = "e) Slab Types"
    ax5.xlabel = "EC steel [kgCO2e/m²]"
    ax5.ylabel = "EC RC-slab [kgCO2e/m²]"

    hull = andrew_hull(orth_biaxial_data.steel_ec, orth_biaxial_data.slab_ec)[1]
    points = [Point2f(orth_biaxial_data.steel_ec[i], orth_biaxial_data.slab_ec[i]) for i in hull]
    poly!(ax5, points, color=(色[:irispurple], 0.05), transparency=true, strokewidth=0.5, strokecolor=色[:irispurple], linestyle=:dash)

    hull = andrew_hull(uniaxial_data.steel_ec, uniaxial_data.slab_ec)[1]
    points = [Point2f(uniaxial_data.steel_ec[i], uniaxial_data.slab_ec[i]) for i in hull]
    poly!(ax5, points, color=(色[:magenta], 0.02), transparency=true, strokewidth=0.5, strokecolor=色[:magenta], linestyle=:dash)

    hull = andrew_hull(isotropic_data.steel_ec, isotropic_data.slab_ec)[1]
    points = [Point2f(isotropic_data.steel_ec[i], isotropic_data.slab_ec[i]) for i in hull]
    poly!(ax5, points, color=(色[:skyblue], 0.05), transparency=true, strokewidth=0.5, strokecolor=色[:skyblue], linestyle=:dash)


    # Add legends
    axislegend(ax1, [sk[:elem_usual], sk[:elem_optimal]], ["Uniform", "Cellular"], position=:rt, orientation=:vertical, labelhalign=:left, framevisible=true, backgroundcolor=(:white, 0), framecolor=(:white, 0), labelsize=sk[:fontsize_small], patchsize=(2, 10), padding=(0, 0, 0, 0))
    axislegend(ax2, [sk[:elem_usual], sk[:elem_optimal]], ["Catalog (W)", "Continuous"], position=:rt, orientation=:vertical, labelhalign=:left, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=sk[:fontsize_small], patchsize=(2, 10), padding=(0, 0, 0, 0))
    axislegend(ax3, [sk[:elem_usual], sk[:elem_optimal]], ["Collinear", "Noncollinear"], position=:rt, orientation=:vertical, labelhalign=:left, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=sk[:fontsize_small], patchsize=(2, 10), padding=(0, 0, 0, 0))
    axislegend(ax4, [sk[:elem_usual], sk[:elem_optimal]], ["0.125m", "0.001m"], position=:rt, orientation=:vertical, labelhalign=:left, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=sk[:fontsize_small], patchsize=(2, 10), padding=(0, 0, 0, 0))
    axislegend(ax5, [sk[:elem_isotropic], sk[:elem_orthogonal], sk[:elem_uniaxial]], ["Isotropic", "Biaxial Orthogonal", "Uniaxial"], position=:rt, orientation=:vertical, labelhalign=:left, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=sk[:fontsize_small], patchsize=(2, 10), padding=(0, 0, 0, 0))

    # Find the business-as-usual (BAU) scenario
    bau_filter = row -> row.name == "r1c2" && row.slab_type == "uniaxial" && row.slab_sizer == "uniform" && row.beam_sizer == "discrete" && row.collinear == true && row.vector_1d_x == 1 && row.vector_1d_y == 0 && row.max_depth == 40
    business_as_usual = filter(bau_filter, df_all)

    bau_steel = business_as_usual.steel_ec
    bau_slab = business_as_usual.concrete_ec + business_as_usual.rebar_ec
    bau_total = bau_steel[1] + bau_slab[1]

    # Add BAU lines to all axes
    for ax in [ax1, ax2, ax3, ax4, ax5]
        vlines!(ax, bau_steel, color=:black, linestyle=:dash, transparency=true, linewidth=1)
        hlines!(ax, bau_slab, color=:black, linestyle=:dash, transparency=true, linewidth=1)
        lines!(ax, [bau_total, 0], [0, bau_total], color=:black, linestyle=:dash, transparency=true, linewidth=1)
        if ax == ax5
            text!(ax, 149, bau_slab[1] + 1, text="BAU slab", color = :black, align = (:right, :bottom), fontsize = sk[:fontsize_small])
            text!(ax, bau_steel[1] + 3, 1, text="BAU steel", color = :black, rotation = pi/2, align = (:left, :center), fontsize = sk[:fontsize_small])
            text!(ax, 3, bau_slab[1]+bau_steel[1]+1, text="BAU total", color = :black, rotation=-pi/4, align = (:left, :center), fontsize = sk[:fontsize_small])    
        end
    end

    # Link axes
    linkxaxes!(ax1, ax2, ax3, ax4, ax5)
    linkyaxes!(ax1, ax2, ax3, ax4, ax5)

    di = DataInspector(fig)

    GC.gc()

    return fig

end
