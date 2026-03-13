function plot_6_depth(df)

    depths = sort(unique(df.max_depth))
    slab_names = unique(df.name)
    fontsize = 11

    slab_sizer = "uniform"
    colors = [色[:skyblue],色[:magenta],色[:skyblue]]
    beam_sizers = ["discrete", "continuous", "discrete"]
    alphas = [1,1,0.25]
    collinear = true

    elem_discrete = MarkerElement(color = 色[:skyblue], marker = :star8)
    elem_continuous = MarkerElement(color = 色[:magenta], marker = :star8)

    fig = Figure(size=(190*4,190*2))
    ax = Axis(fig[1,1], title="Maximum depths", xlabel = "Maximum assembly depth [in]", ylabel = "Steel EC [kgCO2e/m²]", limits = (0,60,0,nothing), topspinevisible=false, rightspinevisible=false, xticks = (0:10:60, ["0", "10", "20", "30", "40", "50", "∞"]), yticklabelsize = fontsize, xticklabelsize = fontsize, xlabelsize = fontsize, ylabelsize = fontsize, titlesize=fontsize)

    for i in 1:lastindex(beam_sizers)

        beam_sizer = beam_sizers[i]
        color = colors[i]
        alpha = alphas[i]

        for name in slab_names

            y = Float64[]
            
            for depth in depths
                filtered_slab = filter(row -> row.name == name && row.slab_sizer == slab_sizer && row.max_depth == depth && row.collinear == collinear && row.beam_sizer == beam_sizer, df)
                if length(filtered_slab.name) == 0; push!(y, NaN); continue; end
                push!(y, filtered_slab.steel_ec[1])
                #push!(y, min(filtered_slab.steel_ec[1], minimum(filter(!isnan, y))))
            end

            lines!(ax, depths, y, color = color, alpha=alpha, transparency=true)
            scatter!(ax, depths, y, color = color, alpha=alpha, transparency=true, marker=:star8)

        end

    end

    axislegend(ax, [elem_discrete, elem_continuous], ["Discrete", "Continuous"], position=:lb, orientation = :vertical, labelhalign = :left, framevisible = true, framecolor=:white, backgroundcolor= :white, labelsize = fontsize)

    GC.gc()

    return fig

end
