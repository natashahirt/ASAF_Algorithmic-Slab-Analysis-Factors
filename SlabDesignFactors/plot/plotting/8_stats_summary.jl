function plot_8_stats_summary(df_all)

    df_all = filter(row -> !isnan(row.total_ec), df_all)

    filter_slab_sizer = row -> row.slab_sizer == "cellular"
    filter_beam_sizer = row -> row.beam_sizer == "continuous"
    filter_collinearity = row -> row.collinear == false
    filter_max_depth = row -> row.max_depth == 40
    filter_slab_type = row -> row.slab_type == "uniaxial"
    filter_type_and_vector = row -> row.slab_type == "uniaxial" && row.vector_1d_x == 1 && row.vector_1d_y == -1
    filter_category = row -> row.category == "nova"

    filter_conventional = row -> row.slab_type == "uniaxial" && row.vector_1d_x == 1 && row.vector_1d_y == 0 && row.beam_sizer == "discrete" && row.slab_sizer == "uniform" && row.collinear == true && row.max_depth == 40
    filter_bau = row -> row.name == "r1c2" && row.slab_type == "uniaxial" && row.slab_sizer == "uniform" && row.beam_sizer == "discrete" && row.collinear == true && row.vector_1d_x == 1 && row.vector_1d_y == 0 && row.max_depth == 40
    df_below_bau = filter(row -> row.total_ec < filter(filter_bau,df_all).total_ec[1], df_all)
    
    df_filtered = filter(filter_bau, df_all)

    round(mean(df_filtered.steel_ec), digits=2)
    round(mean(df_filtered.slab_ec), digits=2)
    round(mean(df_filtered.total_ec), digits=2)
    round(std(df_filtered.steel_ec), digits=2)
    round(std(df_filtered.slab_ec), digits=2)
    round(std(df_filtered.total_ec), digits=2)

    round(mean(df_below_bau.steel_ec), digits=2)
    round(mean(df_below_bau.slab_ec), digits=2)
    round(mean(df_below_bau.total_ec), digits=2)
    round(std(df_below_bau.steel_ec), digits=2)
    round(std(df_below_bau.slab_ec), digits=2)
    round(std(df_below_bau.total_ec), digits=2)

    bau = filter(filter_bau,df_all).total_ec[1] # business as usual
    reduction_conventional = (bau - minimum(filter(filter_conventional, df_all).total_ec)) / bau # reduction for conventional
    reduction_all = (bau - minimum(df_all.total_ec)) / bau # reduction for all 

    ###### ============================================================

    fontsize = 11

    fig = Figure(size=(190*4,190*2))
    master_grid = GridLayout(fig[1,1])
    grids = [GridLayout(master_grid[1,1]), GridLayout(master_grid[1,2])]
    colsize!(master_grid,1,Relative(1/2))
    texts = ["a) Full dataset", "b) Total EC < business-as-usual"]

    for (k, df_plot) in enumerate([df_all, df_below_bau])

        grid = grids[k]
        Label(grid[0, :], text = texts[k], fontsize = fontsize, font = :bold, tellwidth = false)

        ax1 = Axis(grid[1,1], title = "Total", ylabel = "EC [kgCO2e/m²]", xticks = (1:5, ["Slab sizing", "Beam\nsizing", "Beam\ncollinearity", "Assembly\ndepth", "Slab types"]), limits = (nothing,nothing,0,150), titlesize = fontsize, yticklabelsize = fontsize, xticklabelsize = fontsize, xlabelsize = fontsize, ylabelsize = fontsize)
        ax2 = Axis(grid[2,1], title = "Steel", ylabel = "EC [kgCO2e/m²]", xticks = (1:5, ["Slab sizing", "Beam\nsizing", "Beam\ncollinearity", "Assembly\ndepth", "Slab types"]), limits = (nothing,nothing,0,150), titlesize = fontsize, yticklabelsize = fontsize, xticklabelsize = fontsize, xlabelsize = fontsize, ylabelsize = fontsize)
        ax3 = Axis(grid[3,1], title = "Slab", ylabel = "EC [kgCO2e/m²]", xticks = (1:5, ["Slab sizing", "Beam\nsizing", "Beam\ncollinearity", "Assembly\ndepth", "Slab types"]), limits = (nothing,nothing,0,150), titlesize = fontsize, yticklabelsize = fontsize, xticklabelsize = fontsize, xlabelsize = fontsize, ylabelsize = fontsize)

        dodge_slab_sizer = [df_plot[i,:].slab_sizer == "uniform" ? 1 : 2 for i in 1:lastindex(df_plot.name)]
        dodge_beam_sizer = [df_plot[i,:].beam_sizer == "discrete" ? 1 : 2 for i in 1:lastindex(df_plot.name)]
        dodge_collinear = [df_plot[i,:].collinear == true ? 1 : 2 for i in 1:lastindex(df_plot.name)]
        dodge_max_depth = [df_plot[i,:].max_depth == 25 ? 1 : 2 for i in 1:lastindex(df_plot.name)]
        dodge_slab_type = [df_plot[i,:].slab_type == "isotropic" ? 1 : df_plot[i,:].slab_type == "orth_biaxial" ? 2 : 3 for i in 1:lastindex(df_plot.name)]

        for (i,dodge) in enumerate([dodge_slab_sizer, dodge_beam_sizer, dodge_collinear, dodge_max_depth, dodge_slab_type])

            if length(unique(dodge)) == 2
                color_map = map(d -> d == 1 ? 色[:skyblue] : 色[:irispurple], dodge)
            else
                color_map = map(d -> d == 1 ? 色[:skyblue] : d == 2 ? 色[:irispurple] : 色[:magenta], dodge)
            end

            for (j,ax) in enumerate([ax1, ax2, ax3])
                data = [df_plot.total_ec, df_plot.steel_ec, df_plot.slab_ec][j]
                boxplot!(ax, ones(length(df_plot.name)) * i, data, dodge = dodge, color = color_map, mediancolor=:black, whiskerlinewidth=0.5, medianlinewidth=0.5, markersize = 2)
            end

        end

        linkyaxes!(ax1,ax2,ax3)

    end

    sample_isometric_point = MarkerElement(color = 色[:skyblue], marker = :rect, strokecolor = :black)
    sample_orth_point = MarkerElement(color = 色[:irispurple], marker = :rect, strokecolor = :black)
    sample_uniaxial_point = MarkerElement(color = 色[:magenta], marker = :rect, strokecolor = :black)

    Legend(fig[2,1],[sample_isometric_point, sample_orth_point, sample_uniaxial_point],
        ["Isotropic", "Biaxial Orthogonal", "Uniaxial"], 
        orientation=:horizontal,
        tellwidth=false,
        tellheight=true,
        framevisible=false,
        labelsize=fontsize)

    GC.gc()

    return fig

end


"""
    create_slab_summary_table(df)

Creates a summary table comparing slabs with and without the "drawn" attribute.
The table includes mean values for steel, slab, and total EC for each slab name.
"""
function create_slab_summary_table(df)
    df = filter(row -> row.max_depth == 25, df)

    slab_names = ["triple_bay_bau", "triple_bay", "s6-6"]

    println("slab             standard                                   modified for release")
    println("             -------------------------                  ---------------------------")
    println("             steel ec    slab ec    total ec          steel ec    slab ec    total ec")
    println("             ---------   --------   ---------         ---------   --------   ---------")

    for slab_name in slab_names
        df_slab = filter(row -> row.name == slab_name, df)
        df_slab_drawn = filter(row -> row.name == slab_name * "_drawn" in df.name, df)
        
        # Get mean values for standard slab
        mean_steel_ec = mean(df_slab.steel_ec)
        mean_slab_ec = mean(df_slab.slab_ec) 
        mean_total_ec = mean(df_slab.total_ec)

        mean_steel_ec_drawn = mean(df_slab_drawn.steel_ec)
        mean_slab_ec_drawn = mean(df_slab_drawn.slab_ec)
        mean_total_ec_drawn = mean(df_slab_drawn.total_ec)

        # Print formatted row
        println(rpad(slab_name, 15) * lpad(round(mean_steel_ec, digits=1), 9) * "   " * lpad(round(mean_slab_ec, digits=1), 8) * "   " * lpad(round(mean_total_ec, digits=1), 9) * "         " * lpad(round(mean_steel_ec_drawn, digits=1), 9) * "   " * lpad(round(mean_slab_ec_drawn, digits=1), 8) * "   " * lpad(round(mean_total_ec_drawn, digits=1), 9))

    end
end

"""
df = assemble_data("SlabDesignFactors/results/test_results/constructability.csv")

# Get maximum and minimum total embodied carbon
max_ec = maximum(df_combined.steel_ec)
min_ec = minimum(df_combined.steel_ec)
# Get the row with minimum total embodied carbon
min_ec_row = df_combined[argmin(df_combined.steel_ec), :]

bau_slab.total_ec

# Calculate percentage improvement from BAU
pct_improvement = (bau_slab.steel_ec - min_ec) / bau_slab.total_ec * 100
println("\nPercentage improvement from BAU: ", round(pct_improvement, digits=1), "%")"""


"""
    create_filtered_summary_tables(df)

Creates separate summary tables for each design decision using filtering.
Each table includes mean and standard deviation for steel, slab, and total EC.
"""
function create_summary_tables(df)
    slab_filter = row -> row.name == "r1c2" && row.slab_type == "uniaxial" && row.beam_sizer == "discrete" && row.vector_1d_x == 1 && row.vector_1d_y == 0 && row.slab_sizer == "uniform" && row.max_depth == 40 && row.collinear == true && row.slab_min == true
    bau_slab = filter(slab_filter, df)[1,:]
    df = filter(row -> row.total_ec <= bau_slab.total_ec, df)

    # Define the design decisions, including slab types and vector combinations
    design_decisions = [:slab_sizer, :beam_sizer, :collinear, :slab_min, :slab_type, :max_depth, :category]

    for decision in design_decisions
        println("\n\nDesign Decision: ", decision)
        println("        Variation |             Mean / kgCO₂/m²     | Standard Dev / kgCO₂/m² | Slab count")
        println("                 |              -----  -----  ----- | -----  -----  -----  |  -----")

        # Get unique values for the current design decision
        unique_values = unique(df[!, decision])

        for value in unique_values
            # Filter the dataframe for the current design decision value
            df_filtered = filter(row -> row[decision] == value, df)

            # Calculate statistics
            mean_steel = mean(df_filtered.steel_ec)
            mean_slab = mean(df_filtered.slab_ec)
            mean_total = mean(df_filtered.total_ec)
            std_steel = std(df_filtered.steel_ec)
            std_slab = std(df_filtered.slab_ec)
            std_total = std(df_filtered.total_ec)
            slab_count = nrow(df_filtered)

            # Print the results with design decision in left cell
            println(rpad(string(decision) * ": " * string(value), 30) * 
                    lpad(round(mean_steel, digits=2), 6) * " " *
                    lpad(round(mean_slab, digits=2), 6) * " " *
                    lpad(round(mean_total, digits=2), 6) * " " *
                    lpad(round(std_steel, digits=2), 6) * " " *
                    lpad(round(std_slab, digits=2), 6) * " " *
                    lpad(round(std_total, digits=2), 6) * " " *
                    lpad(string(slab_count), 6))
        end
    end

    # Iterate over slab types and vector combinations
    slab_types = ["isotropic", "orth_biaxial", "orth_biaxial", "uniaxial", "uniaxial", "uniaxial", "uniaxial"]
    vector_combinations = [[0, 0], [1, 0], [1, 1], [1, 0], [0, 1], [1, 1], [1, -1]]

    println("\n\nDesign Decision: ", :slab_type)
    println("        Variation |             Mean / kgCO₂/m²     | Standard Dev / kgCO₂/m² | Slab count")
    println("                 |              -----  -----  ----- | -----  -----  -----  |  -----")

    
    for (slab_type, vector) in zip(slab_types, vector_combinations)

        df_filtered = filter(row -> row.slab_type == slab_type && 
                                          row.vector_1d_x == vector[1] && 
                                          row.vector_1d_y == vector[2], df)

         # Calculate statistics
         mean_steel = mean(df_filtered.steel_ec)
         mean_slab = mean(df_filtered.slab_ec)
         mean_total = mean(df_filtered.total_ec)
         std_steel = std(df_filtered.steel_ec)
         std_slab = std(df_filtered.slab_ec)
         std_total = std(df_filtered.total_ec)
         slab_count = nrow(df_filtered)

         # Print the results with design decision in left cell
         println(rpad(string(slab_type) * " | " * string(vector), 30) * 
                 lpad(round(mean_steel, digits=2), 6) * " " *
                 lpad(round(mean_slab, digits=2), 6) * " " *
                 lpad(round(mean_total, digits=2), 6) * " " *
                 lpad(round(std_steel, digits=2), 6) * " " *
                 lpad(round(std_slab, digits=2), 6) * " " *
                 lpad(round(std_total, digits=2), 6) * " " *
                 lpad(string(slab_count), 6))
    end

    # Print statistics for entire database
    println("\n\nEntire Database Statistics")
    println("        Variation |             Mean / kgCO₂/m²     | Standard Dev / kgCO₂/m² | Slab count")
    println("                 |              -----  -----  ----- | -----  -----  -----  |  -----")

    # Calculate statistics for full dataset
    mean_steel = mean(df.steel_ec)
    mean_slab = mean(df.slab_ec) 
    mean_total = mean(df.total_ec)
    std_steel = std(df.steel_ec)
    std_slab = std(df.slab_ec)
    std_total = std(df.total_ec)
    slab_count = nrow(df)

    # Print the results
    println(rpad("All designs", 30) * 
            lpad(round(mean_steel, digits=2), 6) * " " *
            lpad(round(mean_slab, digits=2), 6) * " " *
            lpad(round(mean_total, digits=2), 6) * " " *
            lpad(round(std_steel, digits=2), 6) * " " *
            lpad(round(std_slab, digits=2), 6) * " " *
            lpad(round(std_total, digits=2), 6) * " " *
            lpad(string(slab_count), 6))

    # Print statistics for BAU slab
    println("\n\nBusiness As Usual Slab Statistics")
    println("        Variation |             Mean / kgCO₂/m²     | Standard Dev / kgCO₂/m² | Slab count")
    println("                 |              -----  -----  ----- | -----  -----  -----  |  -----")

    # Print the results for single BAU slab
    println(rpad("BAU slab", 30) * 
            lpad(round(bau_slab.steel_ec, digits=2), 6) * " " *
            lpad(round(bau_slab.slab_ec, digits=2), 6) * " " *
            lpad(round(bau_slab.total_ec, digits=2), 6) * " " *
            lpad("0.00", 6) * " " *  # Standard deviation is 0 for single value
            lpad("0.00", 6) * " " *
            lpad("0.00", 6) * " " *
            lpad("1", 6))

end



# Example usage
# df = DataFrame(...) # Load your combined dataset here
# create_filtered_summary_tables(df)