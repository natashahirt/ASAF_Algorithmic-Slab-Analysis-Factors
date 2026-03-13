"""
    plot_2_megaplot(df_all)

Creates a megaplot figure to visualize slab design factors with various comparisons.
"""
function plot_2_megaplot(df_all)

    fig = Figure(size=(190*4, 190*4.5))

    # Figure settings
    alpha = 0.8
    alpha_fade = 0.5
    transparency = true
    markersize = 5
    markersize_zoom = 7
    fontsize = 11
    smallfontsize = 8

    # Legends
    elem_grey = MarkerElement(color = :lightgrey, marker = :circle)
    elem_selected = MarkerElement(color = 色[:magenta], marker = :circle)
    elem_isotropic = MarkerElement(color = 色[:ceruleanblue], marker = :nova8)
    elem_orthogonal = MarkerElement(color = 色[:irispurple], marker = :cross)
    elem_uniaxial = MarkerElement(color = 色[:magenta], marker = :hline)
    elem_t = MarkerElement(color = 色[:magenta], marker = :circle)
    elem_g = MarkerElement(color = 色[:lilac], marker = :circle)
    elem_s = MarkerElement(color = :navy, marker = :circle)
    elem_concrete = MarkerElement(color = 色[:charcoalgrey], marker = :circle)
    elem_steel = MarkerElement(color = 色[:skyblue], marker = :circle)
    elem_rebar = MarkerElement(color = 色[:magenta], marker = :circle)

    # Grid
    grid = fig[1, 1] = GridLayout()
    scattergrid = grid[1,1] = GridLayout()
    gradient_grid = scattergrid[2,1] = GridLayout()
    bargrid = grid[1,2] = GridLayout()

    rowsize!(scattergrid, 2, Fixed(203))
    rowgap!(scattergrid, 1, Fixed(20))
    colsize!(grid, 1, Relative(5/6))

    # Import CSVs
    df_all = filter(row -> !isnan(row.total_ec), df_all)

    max_steel = maximum(df_all.steel_ec) * 1.25
    max_slab = maximum(df_all.slab_ec) * 1.25

    # Pick out the individuals, novating with business as usual
    bau_filter = row -> row.name == "r1c2" && row.slab_type == "uniaxial" && row.slab_sizer == "uniform" && row.beam_sizer == "discrete" && row.collinear == true && row.vector_1d_x == 1 && row.vector_1d_y == 0 && row.max_depth == 40
    business_as_usual = filter(bau_filter, df_all)

    bau_steel = business_as_usual.steel_ec[1]
    bau_slab = business_as_usual.slab_ec[1]

    color_filter = row -> row.beam_sizer == "discrete" && row.collinear == true && row.slab_sizer == "uniform"
    df_grey = filter(!color_filter, df_all)
    df_color = filter(color_filter, df_all)

    # Separate larger and smaller than
    filter_function = row -> row.steel_ec + row.slab_ec > bau_steel + bau_slab

    greater_than_bau_grey = filter(filter_function, df_grey)
    less_than_bau_grey = filter(!filter_function, df_grey)
    greater_than_bau_color = filter(filter_function, df_color)
    less_than_bau_color = filter(!filter_function, df_color)

    less_than_bau = filter(!filter_function, df_all)
    # Megaplot -- less than
    x_lim_max = 150
    y_lim_max =150
    ax1 = Axis(scattergrid[1,1], xticks = collect(0:10:maximum(df_all.steel_ec)), title = "a) Full Dataset", aspect = DataAspect(), xlabel = "EC steel [kgCO2e/m²]", ylabel = "EC RC-slab [kgCO2e/m²]", limits=(0,x_lim_max,0,y_lim_max), yticklabelsize = fontsize, xticklabelsize = fontsize, xlabelsize = fontsize, ylabelsize = fontsize, titlesize = fontsize)

    x_lim_min = minimum(less_than_bau.steel_ec) * .3
    y_lim_min = minimum(less_than_bau.slab_ec) * .9

    ax2_limits = (x_lim_min,120,y_lim_min,bau_slab*1.1) #(x_lim_min,bau_steel * 1.1,y_lim_min,bau_slab*1.1)
    ax2 = Axis(gradient_grid[1,1], yticks = collect(0:10:bau_slab), aspect=DataAspect(), title = "b) Element density (inset)", xlabel = "EC steel [kgCO2e/m²]", ylabel = "EC RC-slab [kgCO2e/m²]", yticklabelsize = fontsize, xticklabelsize = fontsize, xlabelsize = fontsize, ylabelsize = fontsize, titlesize=fontsize, limits=ax2_limits)
    axislegend(ax1, [elem_t, elem_g, elem_s, elem_grey], ["Topology", "Grid", "Nova", "Non-standard"], position = :rb, orientation = :vertical, labelhalign = :right, framevisible = true, backgroundcolor= :white, framecolor = :white, labelsize=9, patchsize = (2,10), padding=(0,0,0,0))
    axislegend(ax2, [elem_g, elem_s], ["Grid", "Nova"], position = :rb, orientation = :vertical, labelhalign = :right, framevisible = true, backgroundcolor= :white, framecolor = :white, labelsize=9, patchsize = (2,10), padding=(2,2,2,2))

    ax3 = Axis(bargrid[1,1], title = "c) Ranked EC", xlabel = "EC [kgCO2e/m²]", limits=(0,maximum(df_all.total_ec),0,length(df_all.name)), yticklabelsize = fontsize, xticklabelsize = fontsize, xlabelsize = fontsize, ylabelsize = fontsize, titlesize= fontsize)

    # Iso-EC lines
    ec_lines = collect(range(1, maximum(df_all.total_ec) / 10)) .* 10

    for i in 1:lastindex(ec_lines)
        if ec_lines[i] <= bau_steel + bau_slab
            lines!(ax1, [0, ec_lines[i]], [ec_lines[i], 0], color = ec_lines[i], colormap = (Reverse(:grays), 0.25), colorrange = (0, bau_steel + bau_slab), linestyle = :dash, linewidth = 1, inspectable = false)
            lines!(ax2, [0, ec_lines[i]], [ec_lines[i], 0], color = ec_lines[i], colormap = (Reverse(:grays), 0.25), colorrange = (0, bau_steel + bau_slab), linestyle = :dash, linewidth = 1, inspectable = false)
        else
            lines!(ax1, [0, ec_lines[i]], [ec_lines[i], 0], color = ec_lines[i], colormap = (:grays, 0.25), colorrange = (bau_steel + bau_slab, maximum(ec_lines)), linestyle = :dash, linewidth = 1, inspectable = false)
            lines!(ax2, [0, ec_lines[i]], [ec_lines[i], 0], color = ec_lines[i], colormap = (:grays, 0.25), colorrange = (bau_steel + bau_slab, maximum(ec_lines)), linestyle = :dash, linewidth = 1, inspectable = false)
        end
    end
    
    # Plot inset
    limit_points = [Point2f(ax2_limits[1], ax2_limits[3]), Point2f(ax2_limits[2], ax2_limits[3]), Point2f(ax2_limits[2], ax2_limits[4]), Point2f(ax2_limits[1], ax2_limits[4])]
    poly!(ax1, limit_points, color = :white, inspectable = false, strokewidth = 1, strokecolor = :lightgrey, transparency = false)
    text!(ax1, ax2_limits[1]+1, ax2_limits[3]+1, text="inset", color = :lightgrey, align = (:left, :bottom), fontsize = smallfontsize)

    # Plot greater than
    filter_function = row -> row.category == "topology"
    df_topology = filter(filter_function, greater_than_bau_color)
    filter_function = row -> row.category == "grid"
    df_grid = filter(filter_function, greater_than_bau_color)
    filter_function = row -> row.category == "nova"
    df_nova = filter(filter_function, greater_than_bau_color)

    scatter!(ax1, greater_than_bau_grey.steel_ec, greater_than_bau_grey.slab_ec, marker=greater_than_bau_grey.symbol, rotation=greater_than_bau_grey.rotation, color=(:lightgrey,alpha_fade), transparency = false, markersize=markersize, inspector_label = (self, i, p) -> greater_than_bau_grey.rowcol[i])

    scatter!(ax1, df_topology.steel_ec, df_topology.slab_ec, marker=df_topology.symbol, rotation=df_topology.rotation, color=(色[:magenta],alpha_fade), transparency = false, markersize=markersize, inspector_label = (self, i, p) -> df_topology.rowcol[i])
    scatter!(ax1, df_grid.steel_ec, df_grid.slab_ec, marker=df_grid.symbol, rotation=df_grid.rotation, color=(色[:lilac],alpha_fade), transparency = false, markersize=markersize, inspector_label = (self, i, p) -> df_grid.rowcol[i])
    scatter!(ax1, df_nova.steel_ec, df_nova.slab_ec, marker=df_nova.symbol, rotation=df_nova.rotation, color=(:navy,alpha_fade), transparency = false, markersize=markersize, inspector_label = (self, i, p) -> df_nova.rowcol[i])

    # Plot less than
    filter_function = row -> row.category == "topology"
    df_topology = filter(filter_function, df_color)
    filter_function = row -> row.category == "grid"
    df_grid = filter(filter_function, df_color)
    filter_function = row -> row.category == "nova"
    df_nova = filter(filter_function, df_color)

    scatter!(ax1, df_grey.steel_ec, df_grey.slab_ec, marker=df_grey.symbol, rotation=df_grey.rotation, color=(:lightgrey,alpha), transparency = false, markersize=markersize, inspector_label = (self, i, p) -> df_grey.rowcol[i])

    scatter!(ax1, df_topology.steel_ec, df_topology.slab_ec, marker=df_topology.symbol, rotation=df_topology.rotation, color=(色[:magenta],alpha), transparency = transparency, markersize=markersize, inspector_label = (self, i, p) -> df_topology.rowcol[i])
    scatter!(ax1, df_grid.steel_ec, df_grid.slab_ec, marker=df_grid.symbol, rotation=df_grid.rotation, color=(色[:lilac],alpha), transparency = transparency, markersize=markersize, inspector_label = (self, i, p) -> df_grid.rowcol[i])
    scatter!(ax1, df_nova.steel_ec, df_nova.slab_ec, marker=df_nova.symbol, rotation=df_nova.rotation, color=(:navy,alpha), transparency = transparency, markersize=markersize, inspector_label = (self, i, p) -> df_nova.rowcol[i])

    # Plot the same on ax2
    scatter!(ax2, less_than_bau_grey.steel_ec, less_than_bau_grey.slab_ec, marker=less_than_bau_grey.symbol, rotation=less_than_bau_grey.rotation, color=(:lightgrey,alpha), transparency = false, markersize=markersize_zoom, inspectable = false)
    scatter!(ax2, df_topology.steel_ec, df_topology.slab_ec, marker=df_topology.symbol, rotation=df_topology.rotation, color=(:lightgrey,alpha), transparency = false, markersize=markersize_zoom, inspector_label = (self, i, p) -> df_topology.rowcol[i])

    n_elements_grid = [length(df_grid[i,:].ids) for i in 1:lastindex(df_grid.name)]
    n_elements_nova = [length(df_nova[i,:].ids) for i in 1:lastindex(df_nova.name)]

    try
        scatter!(ax2, df_grid.steel_ec, df_grid.slab_ec, marker=df_grid.symbol, rotation=df_grid.rotation, color = n_elements_grid, colormap = Reverse(:acton), colorrange = (0,maximum(n_elements_nova)), transparency = false, markersize=markersize_zoom, inspector_label = (self, i, p) -> df_grid.rowcol[i])
        Colorbar(gradient_grid[1,2], limits = (0, maximum(n_elements_nova)), colormap = Reverse(:acton), ticklabelsize=smallfontsize, labelsize=smallfontsize, label="# grid elements", ticklabelrotation=pi/2)
    catch
        println("Error plotting grid elements")
    end
    try
        scatter!(ax2, df_nova.steel_ec, df_nova.slab_ec, marker=df_nova.symbol, rotation=df_nova.rotation, color = n_elements_nova, colormap = Reverse(:devon), colorrange = (0,maximum(n_elements_nova)), transparency = false, markersize=markersize_zoom, inspector_label = (self, i, p) -> df_nova.rowcol[i])
        Colorbar(gradient_grid[1,3], limits = (0, maximum(n_elements_nova)), colormap = Reverse(:devon), ticklabelsize=smallfontsize, labelsize=smallfontsize, label="# nova elements", ticklabelrotation=pi/2)
    catch
        println("Error plotting nova elements")
    end

    # Plot business as usual
    scatter!(ax1, bau_steel, bau_slab, marker=business_as_usual.symbol, transparency = transparency, markersize=markersize, rotation=business_as_usual.rotation, color=:black, inspector_label = (self, i, p) -> business_as_usual.rowcol[i])

    vlines!(ax1, bau_steel, color = :black, linestyle = :dash, linewidth = 1, inspectable = false)
    hlines!(ax1, bau_slab, color = :black, linestyle = :dash, linewidth = 1, inspectable = false)
    lines!(ax1, [bau_slab + bau_steel, 0], [0, bau_slab + bau_steel], color = :black, linestyle = :dash, linewidth = 1, inspectable = false)
    text!(ax1, 149, bau_slab + 1, text="BAU slab", color = :black, align = (:right, :bottom), fontsize = smallfontsize)
    text!(ax1, bau_steel + 2, 1, text="BAU steel", color = :black, rotation = pi/2, align = (:left, :center), fontsize = smallfontsize)
    text!(ax1, bau_slab+bau_steel, 3, text="BAU total", color = :black, rotation=-pi/4, align = (:right, :center), fontsize = smallfontsize)

    vlines!(ax2, bau_steel, color = :black, linestyle = :dash, linewidth = 1, inspectable = false)
    hlines!(ax2, bau_slab, color = :black, linestyle = :dash, linewidth = 1, inspectable = false)
    lines!(ax2, [bau_slab + bau_steel, 0], [0, bau_slab + bau_steel], color  = :black, linestyle = :dash, linewidth = 1, inspectable = false)

    # Barplot
    sort!(df_all, order(:total_ec))

    df_less_than = copy(df_all)

    for i in 1:lastindex(df_less_than.name)
        row = df_less_than[i,:]
        if row.total_ec > bau_steel + bau_slab
            row.steel_ec = 0 
            row.concrete_ec = 0
            row.rebar_ec = 0
            row.slab_ec = 0
            row.total_ec = 0
        end
    end
    
    dfs = [df_all, df_less_than]
    alphas = [0.1, 1]

    for i in 1:lastindex(dfs)
        df = dfs[i]
        a = alphas[i]

        positions = 1:lastindex(df.area)
        stack = repeat(positions, inner=3)
        categories = repeat(positions, outer=3)
        colours = repeat([1, 2, 3], inner=lastindex(df.area))

        df_concrete = df.concrete_ec
        df_steel = df.steel_ec
        df_rebar = df.rebar_ec

        data_barplot = vcat(df_concrete, df_steel, df_rebar)

        barplot!(ax3, categories, data_barplot, stack=stack, gap=0, color=colours, direction=:x, colormap = [(色[:charcoalgrey], a), (色[:skyblue], a), (色[:magenta], a)])
    end

    less_than_bau_elements = length(filter(row -> row.total_ec < bau_steel + bau_slab, df_all).name)
    hlines!(ax3, [less_than_bau_elements], color = :black, linestyle = :dash, linewidth = 1, inspectable = false)
    text!(ax3, maximum(df_all.total_ec)-15, less_than_bau_elements-30, text="< BAU", color = :black, align = (:right, :top), fontsize = smallfontsize)
    less_than_axis_limits = length(filter(row -> row.total_ec < 150 + 150, df_all).name)
    hlines!(ax3, [less_than_axis_limits], color = :black, linestyle = :dash, linewidth = 1, inspectable = false)
    text!(ax3, maximum(df_all.total_ec)-15, less_than_axis_limits-30, text="< axis limits", color = :black, align = (:right, :top), fontsize = smallfontsize)
    
    axislegend(ax3, [elem_rebar, elem_steel, elem_concrete], ["Rebar", "Steel", "Concrete"], position = :rb, orientation = :vertical, labelhalign = :right, framevisible = true, backgroundcolor= :white, framecolor= :white, labelsize=smallfontsize, patchsize = (2,10), padding=(0,0,0,0))

    di = DataInspector(fig)

    GC.gc()

    return fig

end