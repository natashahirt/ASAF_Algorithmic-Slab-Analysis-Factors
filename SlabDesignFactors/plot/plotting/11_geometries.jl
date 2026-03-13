function plot_11_geometries(df_all; category="topology")

    category_map = Dict("topology" => "t", "grid" => "g", "nova" => "n")

    df = filter(row -> row.category == category, df_all)
    unique_x = unique(df.row)
    unique_y = unique(df.col)

    ratio = 4/length(unique_y)

    path = "Geometries/$category/"

    size_x = 190*length(unique_y)*ratio
    size_y = 190*length(unique_x)*ratio

    scale = size_x/size_y * 2
    
    fig = Figure(size=(size_y * scale, size_x * scale));

    for i in 1:lastindex(unique_x)
        for j in 1:lastindex(unique_y)
            df_filtered = filter(row -> row.row == unique_x[i] && row.col == unique_y[j], df)
            name = df_filtered.name[1]
            path_name = path * name * ".json"
            geometry_dict = JSON.parse(JSON.parse(replace(read(path_name, String), "\\n" => ""), dicttype=Dict));
            geometry = generate_from_json(geometry_dict, plot=false, drawn=false);
            
            ax = Axis(fig[j, i], aspect=DataAspect(), title=replace(string(df_filtered.rowcol[1]), "$category" => "$(category_map[category])"), titlesize=15, xticklabelsize=15, yticklabelsize=15)
            plot_elements!(ax, geometry.elements, linewidth=0.5)
            hidespines!(ax)
            hidedecorations!(ax, ticklabels=false, ticks=false)
        end
    end

    return fig

end