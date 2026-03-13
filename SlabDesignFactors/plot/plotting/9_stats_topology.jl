function plot_9_stats_topology(df_all)

    df_all = filter(row -> row.category == "topology", df_all)

    filter_bau = row -> row.name == "r1c2" && row.slab_type == "uniaxial" && row.slab_sizer == "uniform" && row.beam_sizer == "discrete" && row.collinear == true && row.vector_1d_x == 1 && row.vector_1d_y == 0 && row.max_depth == 40 && row.slab_min == true
    df_filtered = filter(filter_bau, df_all)
    bau_total_ec = df_filtered.total_ec[1]
    df_below_bau = filter(row -> row.total_ec < bau_total_ec, df_all)

    colors = [色[:magenta], 色[:skyblue], 色[:charcoalgrey]]
    fontsize = 11
    smallfontsize = 7

    fig = Figure(size=(190*4,190*2))
    master_grid = GridLayout(fig[1,1])
    grids = [GridLayout(master_grid[1,1]), GridLayout(master_grid[1,2])]
    colsize!(master_grid,1,Relative(1/2))
    texts = ["a) Full dataset", "b) Total EC < business-as-usual"]
    names = unique(df_all.rowcol)
    labels = replace.(names, "topology" => "")

    for (k, df_considered) in enumerate([df_all, df_below_bau])

        grid = grids[k]
        Label(grid[0, :], text = texts[k], fontsize = fontsize, font = :bold, tellwidth = false)
        df_plot = df_considered
        
        ax1 = Axis(grid[1,1], title = "Total", ylabel = "EC [kgCO2e/m²]", xticks = (1:lastindex(names), labels), limits = (nothing,nothing,0,200), titlesize = fontsize, yticklabelsize = fontsize, xticklabelsize = smallfontsize, xlabelsize = fontsize, ylabelsize = fontsize, xticklabelrotation=pi/2)
        ax2 = Axis(grid[2,1], title = "Steel", ylabel = "EC [kgCO2e/m²]", xticks = (1:lastindex(names), labels), limits = (nothing,nothing,0,200), titlesize = fontsize, yticklabelsize = fontsize, xticklabelsize = smallfontsize, xlabelsize = fontsize, ylabelsize = fontsize, xticklabelrotation=pi/2)
        ax3 = Axis(grid[3,1], title = "Slab", ylabel = "EC [kgCO2e/m²]", xticks = (1:lastindex(names), labels), limits = (nothing,nothing,0,200), titlesize = fontsize, yticklabelsize = fontsize, xticklabelsize = smallfontsize, xlabelsize = fontsize, ylabelsize = fontsize, xticklabelrotation=pi/2)

        # Add horizontal line at business-as-usual total EC level
        hlines!(ax1, bau_total_ec, color=:black, linestyle=:dash, linewidth=1)
        hlines!(ax2, bau_total_ec, color=:black, linestyle=:dash, linewidth=1) 
        hlines!(ax3, bau_total_ec, color=:black, linestyle=:dash, linewidth=1)
        
        for (i, name) in enumerate(names)
            
            df_data = filter(row -> row.rowcol == name, df_plot)

            for (j,ax) in enumerate([ax1, ax2, ax3])
                data = [df_data.total_ec, df_data.steel_ec, df_data.slab_ec][j]
                boxplot!(ax, ones(length(df_data.name)) * i, data, color = colors[j], mediancolor=:black, whiskerlinewidth=1, medianlinewidth=1, markersize = 5)
            end

        end

        linkyaxes!(ax1,ax2,ax3)

    end

    GC.gc()

    return fig

end