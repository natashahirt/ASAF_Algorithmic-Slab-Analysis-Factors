function plot_4_surface(df; category::String="grid")

    # * FIGURE
    alpha = 0.8
    alpha_fade = 0.5
    transparency = true
    markersize = 5
    markersize_zoom = 7
    fontsize = 11
    smallfontsize = 8

    maximum_ec = maximum(filter(x -> !isnan(x), df.total_ec)) * 1.05
    minimum_ec = minimum(filter(x -> !isnan(x), df.total_ec)) * 0.95

    if category == "grid" || category == "nova"
        df_master = filter(row -> row.category == category, df)
    else
        @assert false "Category is required for this plot. Pick one of 'grid' or 'nova'."
    end

    unique_x = unique(df_master.row)
    unique_y = unique(df_master.col)

    azimuth = π/4 #1.275π - π
    
    filter_function = row -> row.slab_sizer == "uniform" && row.collinear == true && row.beam_sizer == "discrete" && row.max_depth == 40 && row.slab_type == "isotropic"
    df = filter(filter_function, df_master)
    x = unique_x
    y = unique_y

    df = add_missing_rows(df, unique_x, unique_y, maximum_ec)

    z = zeros(length(unique_x), length(unique_y))

    for i in 1:lastindex(df.name)
        row = df[i, :]
        z[row.col, row.row] = row.total_ec
    end

    # fine grid
    interp = cubic_spline_interpolation((range(extrema(x)...), range(extrema(y)...)), z, extrapolation_bc=Line()) # interpolation object
    resolution = 50 # multiplication factor

    x2 = range(extrema(x)..., length=length(x)*resolution)
    y2 = range(extrema(y)..., length=length(x)*resolution)
    z2 = [interp(x, y) for x in x2, y in y2]

    # fine grid extrapolated slightly
    x2_extrap = range(minimum(x)-0.5, maximum(x)+0.5, length=length(x)*resolution)
    y2_extrap = range(minimum(y)-0.5, maximum(y)+0.5, length=length(y)*resolution)
    z2_extrap = [interp(x, y) for x in x2_extrap, y in y2_extrap]

    fig = Figure(size=(190*4, 190*3.8))
    grid = GridLayout(fig[1, 1], height=190*3.7)
    grid_top = GridLayout(grid[1, 1])
    grid_surface = GridLayout(grid_top[1, 1])
    grid_sections = GridLayout(grid_top[2, 1])

    rowsize!(grid, 1, Relative(2/3))
    rowsize!(grid_top, 1, Relative(2/3))

    ax1 = Axis3(grid_surface[1, 1], title="a) Design space 3D", xlabel="x", ylabel="y", zlabel="EC [kgCO2e/m²]", zlabeloffset=40, azimuth=azimuth, xticks=x, yticks=y, aspect=(1, 1, 1), yticklabelsize=fontsize, xticklabelsize=fontsize, zticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, zlabelsize=fontsize, titlesize=fontsize, limits=(0, maximum(x), 0, maximum(y), 0, maximum_ec))
    ax2 = Axis(grid_surface[1, 2], title="b) Design space projected", xlabel="x", ylabel="y", aspect=DataAspect(), xticks=x, yticks=y, yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, titlesize=fontsize)

    n_contours = 20

    # 3D surface
    surface!(ax1, x2, y2, z2, colormap=:dense, transparency=true, alpha=alpha, colorrange=[minimum_ec, maximum_ec], highclip=:lightgrey, rasterize=true)
    scatter!(ax1, x, y, z, marker=df.symbol, rotation=df.rotation, color=:black, markersize=markersize, inspector_label=(self, i, p) -> df.rowcol[i])
    contour3d!(ax1, x2, y2, z2, color=:white, linestyle=:dash, alpha=0.5, transparency=true, levels=n_contours, linewidth=0.5)

    # 2D heatmap
    heatmap!(ax2, x2_extrap, y2_extrap, z2_extrap, colormap=:dense, colorrange=[minimum_ec, maximum_ec], highclip=:lightgrey)
    contour!(ax2, x2_extrap, y2_extrap, z2_extrap, color=:white, linestyle=:dash, levels=n_contours, linewidth=0.5)
    scatter!(ax2, df.row, df.col, marker=df.symbol, rotation=df.rotation, color=:black, markersize=markersize, inspector_label=(self, i, p) -> df.rowcol[i])

    ax2.yreversed = true

    Colorbar(grid_surface[1, 3], limits=(minimum(z), maximum(z)), colormap=:dense, tellheight=false, label="EC [kgCO2e/m²]", labelsize=fontsize)

    # LINE SECTIONS

    filter_function = row -> row.slab_sizer == "uniform" && row.collinear == false && row.beam_sizer == "discrete" && row.max_depth == 25 && row.slab_type == "isotropic"
    df = filter(filter_function, df_master)

    x = unique_x
    y = unique_y

    df = add_missing_rows(df, unique_x, unique_y, maximum_ec)

    z_steel = zeros(length(unique_x), length(unique_y))
    z_concrete = zeros(length(unique_x), length(unique_y))
    z_rebar = zeros(length(unique_x), length(unique_y))

    for i in 1:lastindex(df.name)
        row = df[i, :]
        z_steel[row.col, row.row] = row.steel_ec
        z_concrete[row.col, row.row] = row.concrete_ec
        z_rebar[row.col, row.row] = row.rebar_ec
    end

    max_val = maximum([z_steel..., z_concrete..., z_rebar...])
    yticks = collect(0:20:max_val+10)[1:end-1]

    # sample the sections
    section_points_a = [[1, 1], [2, 2], [3, 3], [4, 4], [5, 5], [6, 6]]
    section_points_b = [[1, 6], [2, 5], [3, 4], [4, 3], [5, 2], [6, 1]]

    # Draw lines through section points on ax2
    lines!(ax2, [p[1] for p in section_points_a], [p[2] for p in section_points_a], color=:black, linestyle=:dash, linewidth=1.5, label="Section A")
    lines!(ax2, [p[1] for p in section_points_b], [p[2] for p in section_points_b], color=:black, linestyle=:dash, linewidth=1.5, label="Section B") 
    text!(ax2, section_points_a[end][1], section_points_a[end][2], text="Section A", color=:black, fontsize=fontsize, rotation=-pi/4, align=(:right, :bottom))
    text!(ax2, section_points_b[end][1], section_points_b[end][2], text="Section B", color=:black, fontsize=fontsize, rotation=pi/4, align=(:right, :bottom))

    # Figure
    max_not_nan = maximum(filter(x -> !isnan(x) && x < maximum_ec, df.total_ec))
    yticks_50 = collect(0:50:max_not_nan)
    kwargs = (topspinevisible=false, rightspinevisible=false, yticks=yticks_50, yminorgridvisible=true, yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, titlesize=fontsize)
    ax_a = Axis(grid_sections[1, 1], title="Section A", ylabel="EC [kgCO2e/m²]", limits=(0.5, lastindex(section_points_a)+0.5, 0, max_not_nan), xticks=(collect(1:lastindex(section_points_a)), string.(section_points_a)); kwargs...)
    ax_b = Axis(grid_sections[1, 2], title="Section B", limits=(0.5, lastindex(section_points_b)+0.5, 0, max_not_nan), xticks=(collect(1:lastindex(section_points_b)), string.(section_points_b)); kwargs...)

    x_a = collect(1:length(section_points_a))
    x_b = collect(1:length(section_points_b))

    colors = [色[:skyblue], 色[:charcoalgrey], 色[:magenta]]
    resolution = 50

    for (i, z) in enumerate([z_steel, z_concrete, z_rebar])

        z_a = [z[coord...] for coord in section_points_a]
        z_b = [z[coord...] for coord in reverse(section_points_b)]

        interp_a = cubic_spline_interpolation(range(extrema(x_a)...), z_a, extrapolation_bc=Line())
        interp_b = cubic_spline_interpolation(range(extrema(x_b)...), z_b, extrapolation_bc=Line())

        x_a2 = range(extrema(x_a)..., length=length(x_a)*resolution)
        z_a2 = [interp_a(x) for x in x_a2]

        x_b2 = range(extrema(x_b)..., length=length(x_b)*resolution)
        z_b2 = [interp_b(x) for x in x_b2]

        scatter!(ax_a, x_a, z_a, color=colors[i])
        lines!(ax_a, x_a2, z_a2, color=colors[i])
        scatter!(ax_b, x_b, z_b, color=colors[i])
        lines!(ax_b, x_b2, z_b2, color=colors[i])

    end

    linkyaxes!(ax_a, ax_b)

    elem_concrete = MarkerElement(color=色[:charcoalgrey], marker=:circle)
    elem_steel = MarkerElement(color=色[:skyblue], marker=:circle)
    elem_rebar = MarkerElement(color=色[:magenta], marker=:circle)

    Legend(grid_sections[2, :], [elem_rebar, elem_steel, elem_concrete], ["Rebar", "Steel", "Concrete"], orientation=:horizontal, labelhalign=:right, framevisible=true, backgroundcolor=:white, framecolor=:white, labelsize=fontsize)

    display(fig)

    # HEATMAP SMALL MULTIPLES

    slab_types = ["isotropic", "orth_biaxial", "orth_biaxial", "uniaxial", "uniaxial", "uniaxial", "uniaxial"]
    vector_1ds = [[0., 0.], [1., 0.], [1., 1.], [1., 0.], [0., 1.], [1., 1.], [1., -1.]]
    slab_sizers = ["uniform", "cellular"]

    slab_types_dict = Dict("isotropic" => "Isotropic", "orth_biaxial" => "Orth-Biaxial", "uniaxial" => "Uniaxial")
    slab_sizer_dict = Dict("uniform" => "Uniform", "cellular" => "Cellular")

    grid_smallmul = GridLayout(grid[2, 1], height=160)

    for i in 1:lastindex(slab_sizers)

        slab_sizer = slab_sizers[i]

        for j in 1:lastindex(slab_types)

            vector_1d = vector_1ds[j]
            slab_type = slab_types[j]
            ax = Axis(grid_smallmul[i, j], aspect=1, xticks=x, yticks=y, yticklabelsize=fontsize, xticklabelsize=fontsize, xlabelsize=fontsize, ylabelsize=fontsize, titlesize=fontsize)

            if i == 2
                ax.xlabel = "$(slab_types_dict[slab_type]) \n$(vector_1d)"
            else
                ax.xticklabelsvisible = false
                ax.xticksvisible = false
            end
            if j == 1
                ax.ylabel = slab_sizer_dict[slab_sizer]
                ax.ylabelpadding = 14
            else
                ax.yticklabelsvisible = false
                ax.yticksvisible = false
            end        

            filter_function = row -> row.slab_sizer == slab_sizer && row.slab_type == slab_type && row.vector_1d_x == vector_1d[1] && row.vector_1d_y == vector_1d[2] && row.collinear == true && row.beam_sizer == "discrete" && row.max_depth == 40
            df_smallmul = filter(filter_function, df_master)

            df_smallmul = add_missing_rows(df_smallmul, unique_x, unique_y, maximum_ec)

            make_heatmap(ax, df_smallmul, minimum_ec, maximum_ec, markersize)

        end

    end

    Colorbar(grid_smallmul[:, length(slab_types)+1], limits=(minimum_ec, maximum_ec), colormap=:dense, tellheight=true, label="EC [kgCO2e/m²]", labelsize=fontsize)

    GC.gc()

    return fig

end

function make_heatmap(ax, df_smallmul, minimum_ec, maximum_ec, markersize)

    x = unique(df_smallmul.row)
    y = unique(df_smallmul.col)

    x_incomplete = df_smallmul.row
    y_incomplete = df_smallmul.col
    z_incomplete = df_smallmul.total_ec

    x_complete = repeat(x, inner=lastindex(y))
    y_complete = repeat(y, outer=lastindex(x))
    z_complete = Float64[]

    for (x, y) in collect(zip(x_complete, y_complete))
        k = findfirst(((xi, yi),) -> xi == x && yi == y, collect(zip(x_incomplete, y_incomplete)))

        if k === nothing
            push!(z_complete, maximum_ec + 20)  # Add NaN if the coordinate is missing
        else
            push!(z_complete, z_incomplete[k])  # Add the corresponding z value
        end
    end

    z = reshape(z_complete, (length(x), length(y)))

    # fine grid
    interp = linear_interpolation((range(extrema(x)...), range(extrema(y)...)), z, extrapolation_bc=Line()) # interpolation object
    resolution = 10 # multiplication factor

    # fine grid extrapolated slightly
    x2_extrap = range(minimum(x)-0.5, maximum(x)+0.5, length=length(x)*resolution)
    y2_extrap = range(minimum(y)-0.5, maximum(y)+0.5, length=length(y)*resolution)
    z2_extrap = [interp(x, y) for x in x2_extrap, y in y2_extrap]

    for z in z2_extrap
        if z > maximum_ec
            z = NaN
        end
    end

    n_contours = 10

    # 2D heatmap
    heatmap!(ax, x2_extrap, y2_extrap, z2_extrap, colormap=:dense, colorrange=[minimum_ec, maximum_ec], highclip=:lightgrey)
    contour!(ax, x2_extrap, y2_extrap, z2_extrap, color=:white, linestyle=:dash, levels=minimum_ec:(maximum_ec - minimum_ec)/n_contours:maximum_ec, linewidth=0.3)
    scatter!(ax, x_complete, y_complete, marker=[df_smallmul.symbol[1] for i in 1:lastindex(x_complete)], rotation=[df_smallmul.rotation[1] for i in 1:lastindex(x_complete)], color=:black, inspector_label=(self, i, p) -> df_smallmul.rowcol[i], markersize=markersize)

    ax.yreversed = true
end

function add_missing_rows(df, unique_x, unique_y, maximum_ec)

    for r in unique_x
        for c in unique_y
            if length(filter(row -> row.row == r && row.col == c, df).name) == 0
                row_dict = Dict(name => df[1,name] for name in names(df))
                row_dict["row"] = r
                row_dict["col"] = c 
                row_dict["total_ec"] = maximum_ec + 20
                row_dict["steel_ec"] = maximum_ec + 20
                row_dict["concrete_ec"] = maximum_ec + 20
                row_dict["rebar_ec"] = maximum_ec + 20
                row_dict["rowcol"] = "$r,$c"
                push!(df, row_dict)
            end
        end
    end

    sort!(df, [:row, :col])
    return df

end