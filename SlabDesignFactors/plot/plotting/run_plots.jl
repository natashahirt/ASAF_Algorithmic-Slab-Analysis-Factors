include("../../scripts/_scripts.jl")
include("_plotting.jl")

CairoMakie.activate!()

deflection = "no"
slabmin = "no"

df_all = assemble_data("SlabDesignFactors/results/processed_$(deflection)deflection_$(slabmin)slabmin/")
df_depths = assemble_data("SlabDesignFactors/results/processed_$(deflection)deflection_$(slabmin)slabmin/max_depths_sequential.csv")
df_unfixed = assemble_data("SlabDesignFactors/results/processed_$(deflection)deflection_$(slabmin)slabmin/fix_params.csv")

# get massive multiplot
df_minslabno = assemble_data("SlabDesignFactors/results/processed_$(deflection)deflection_noslabmin/")
df_minslabyes = assemble_data("SlabDesignFactors/results/processed_$(deflection)deflection_yesslabmin/")

df_combined = assemble_data([df_minslabno, df_minslabyes], "slab_min", [false, true])
slab_filter = row -> row.name == "r1c2" && row.slab_type == "uniaxial" && row.beam_sizer == "discrete" && row.vector_1d_x == 1 && row.vector_1d_y == 0 && row.slab_sizer == "uniform" && row.max_depth == 40 && row.collinear == true && row.slab_min == true
bau_slab = filter(slab_filter, df_combined)[1,:]
bau_slab.total_ec

min_ec_row = df_combined[argmin(df_combined.total_ec), :]
min_ec_row.slab_ec
println(min_ec_row)

(96.2-85.54)/96.2

manufacturable_filter = row -> row.beam_sizer == "discrete" && row.slab_sizer == "uniform" && row.collinear == true && row.slab_min == true
manufacturable_slabs = filter(manufacturable_filter, df_combined)
min_ec_row = manufacturable_slabs[argmin(manufacturable_slabs.total_ec), :]
min_ec_row.total_ec
println(min_ec_row)


manufacturable_slabs = filter(manufacturable_filter, df_combined)
min_ec = minimum(manufacturable_slabs.total_ec)
min_ec_row = df_combined[argmin(df_combined.steel_ec), :]

pct_improvement = (bau_slab.total_ec - min_ec) / bau_slab.total_ec * 100

# Count how many slabs have lower embodied carbon than business-as-usual
n_better = count(row -> row.total_ec <= bau_slab.total_ec, eachrow(df_combined))
println("Number of slabs with lower embodied carbon than BAU: $n_better out of $(nrow(df_combined)) total slabs")

save_path = "SlabDesignFactors/plot/figures/$(deflection)deflection/"

fig = plot_1_multiplot(df_combined)
save(save_path * "1_multiplot.pdf", fig)

fig = plot_2_megaplot(df_combined)
save(save_path * "2_megaplot.png", fig)

fig = plot_3_topology(df_combined, category="topology")
save(save_path * "3_topology.pdf", fig)

fig = plot_4_surface(df_all, category="grid")
save(save_path * "4_surface_grid_nominslab_40.pdf", fig)
fig = plot_4_surface(df_all, category="nova")
save(save_path * "4_surface_nova_nominslab_40.pdf", fig)

fig = plot_5_beam_sizes(df_combined, category="topology")
save(save_path * "5_beam_sizes_topology.pdf", fig)
fig = plot_5_beam_sizes(df_combined, category="grid")
save(save_path * "5_beam_sizes_grid.pdf", fig)
fig = plot_5_beam_sizes(df_combined, category="nova")
save(save_path * "5_beam_sizes_nova.pdf", fig)
fig = plot_5_beam_sizes(df_combined)
save(save_path * "5_beam_sizes_all.pdf", fig)

fig = plot_5_beam_sizes_topology(df_combined)
save(save_path * "5_beam_sizes_layouts.pdf", fig)

fig = plot_6_depth(df_depths)
save(save_path * "6_depth_sequential.pdf", fig)

fig = plot_7_fix_params(df_all, df_unfixed)
save(save_path * "7_fix_params.pdf", fig)

fig = plot_8_stats_summary(df_combined)
save(save_path * "8_stats_summary.pdf", fig)

print_8_stats_summary(df_combined)

fig = plot_9_stats_topology(df_combined)
save(save_path * "9_stats_topology.pdf", fig)

fig = plot_10_subplots(df_combined, subplot=:slab_type)
save(save_path * "10_subplots_slab_type.pdf", fig)
fig = plot_10_subplots(df_combined, subplot=:slab_sizer)
save(save_path * "10_subplots_slab_sizer.pdf", fig)
fig = plot_10_subplots(df_combined, subplot=:beam_sizer)
save(save_path * "10_subplots_beam_sizer.pdf", fig)
fig = plot_10_subplots(df_combined, subplot=:collinearity)
save(save_path * "10_subplots_collinearity.pdf", fig)
fig = plot_10_subplots(df_combined, subplot=:max_depth)
save(save_path * "10_subplots_max_depth.pdf", fig)
fig = plot_10_subplots(df_combined, subplot=:slab_min)
save(save_path * "10_subplots_slab_min.pdf", fig)

fig = plot_11_geometries(df_all, category="topology")
save(save_path * "11_geometries_topology.pdf", fig)

# Plot individual slabs

path = "Geometries/grid/x5y5.json"
name = basename(splitext(path)[1])    # Name for the plot
slab_filter = row -> row.name == name && row.slab_type == "isotropic" && row.beam_sizer == "discrete" && row.vector_1d_x == 0 && row.vector_1d_y == 0 && row.slab_sizer == "uniform" && row.max_depth == 40

test_result = filter(slab_filter, df_all)

for row in eachrow(test_result)
    println("Collinear: $(row.collinear), Total EC: $(row.total_ec), Steel EC: $(row.steel_ec)")
end

#df_all_sorted = sort(df_all, :total_ec)
#test_result = slab_params = beam_sizing_params = nothing

try

    #df_slab = filter(slab_filter, df_all)
    test_result = filter(slab_filter, df_all)
    if !(test_result isa DataFrameRow)
        test_result = test_result[1, :] # Convert DataFrame to DataFrameRow by taking first row
    end

    path = "Geometries/$(test_result.category)/$(test_result.name).json"

    # Parse geometry from JSON
    geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict));
    sections = parse_sections(String.(test_result.sections))
    material = material_dict[:in]
    sections = [toASAPframe_W(W_imperial(section), material.E, material.G; convert=false) for section in sections]
    geometry = generate_from_json(geometry_dict, plot=false, drawn=false, sections=sections);

    # Analyze the slab to get dimensions
    slab_params = SlabAnalysisParams(
        geometry, 
        slab_name=test_result.name,
        slab_type=Symbol(test_result.slab_type),
        vector_1d=[test_result.vector_1d_x, test_result.vector_1d_y], 
        slab_sizer=Symbol(test_result.slab_sizer),
        spacing=.1, 
        plot_analysis=true,
        fix_param=true, 
        slab_units=:m,
    );

    if !isnothing(test_result)
        fig = plot_slab(slab_params, test_result, text=false)
    else
        println("No test result found for $name")
        fig = plot_slab(slab_params, beam_sizing_params, text=false)
    end

    save(save_path * "12_$(test_result.name)_$(test_result.collinear).pdf", fig)

catch

    # Parse geometry from JSON
    geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict));
    geometry = generate_from_json(geometry_dict, plot=false, drawn=false);

    # Analyze the slab to get dimensions
    slab_params = SlabAnalysisParams(
        geometry, 
        slab_name=test_result.name,
        slab_type=Symbol(test_result.slab_type),
        vector_1d=[test_result.vector_1d_x, test_result.vector_1d_y], 
        slab_sizer=Symbol(test_result.slab_sizer),
        spacing=.1, 
        plot_analysis=true,
        fix_param=true, 
        slab_units=:m,
    );

    # Sizing parameters
    beam_sizing_params = SlabSizingParams(
        live_load=psf_to_ksi(50), # ksi
        superimposed_dead_load=psf_to_ksi(15), # ksi
        live_factor=1.6, # -
        dead_factor=1.2, # -
        beam_sizer=Symbol(test_result.beam_sizer),
        max_depth=Symbol(test_result.max_depth), # in
        beam_units=:in, # in, etc.
        serviceability_lim=360,
        collinear=Bool(test_result.collinear),
        minimum_continuous=true
    );

    if !isnothing(test_result)
        fig = plot_slab(slab_params, test_result, text=false)
    else
        println("No test result found for $name")
        fig = plot_slab(slab_params, beam_sizing_params, text=false)
    end

    save(save_path * "12_$(test_result.name)_$(test_result.collinear).pdf", fig)

end;


# ==============================
slab_filter = row ->    row.name == "r1c2" && 
                        row.slab_type == "uniaxial" && 
                        row.beam_sizer == "discrete" && 
                        row.vector_1d_x == 1 && 
                        row.vector_1d_y == 0 && 
                        row.slab_sizer == "uniform" && 
                        row.max_depth == 40 && 
                        row.collinear == true

df_slab = filter(slab_filter, df_all)
test_result = df_slab[1, :] # Get first row, can change index as needed

local_deflections, global_deflections = get_max_deflection(test_result);

println("Maximum local deflection: $(maximum(local_deflections)) in")
println("Maximum global deflection: $(maximum(global_deflections)) in")

local_deflections, global_deflections = get_max_deflection(df_combined)
# Filter out NaN values from deflections
local_deflections = filter(!isnan, local_deflections)
global_deflections = filter(!isnan, global_deflections)

println("Maximum local deflection: $(mean(local_deflections)) in")
println("Maximum global deflection: $(mean(global_deflections)) in")

clean_csv_data("SlabDesignFactors/results/processed_nodeflection_noslabmin/max_depths.csv")
clean_csv_data("SlabDesignFactors/results/processed_nodeflection_noslabmin/fix_params.csv")


# ==============================
save_path = "SlabDesignFactors/plot/figures/constrained_inventory/"

df_inventory = assemble_data("SlabDesignFactors/results/constrained_inventory_resized/")
df_max_depth = assemble_data("SlabDesignFactors/results/constrained_inventory_max_depth/")

fig = plot_mass_against_unique_sections(df_inventory)
save(save_path * "4_mass_unique_sections.pdf", fig)

fig = plot_mass_against_max_depth(df_max_depth)
save(save_path * "4_mass_max_depth.pdf", fig)

fig = plot_section_distribution(df_inventory, name="warehouse", max_height=false)

for name in ["office", "school", "warehouse"]
    fig = plot_section_distribution(df_inventory, name=name, max_height=false)
    save(save_path * "4_sections_$(name).pdf", fig)
end

df_filtered = filter(row -> row.max_depth == 25 && row.name == "warehouse", df_max_depth)

# Extract depth and weight from section names using regex
depths = Float64[]
weights = Float64[]

# Convert comma-separated string of sections into vector
sections = split(df_filtered.sections[1], ",")
# Remove any whitespace and brackets
sections = strip.(sections)
sections = replace.(sections, r"[\"\\[\\]]" => "") # Remove quotes and brackets from section names
sections = replace.(sections, "Any" => "") # Remove Any prefix if present
sections = replace.(sections, "\"" => "") # Remove Any prefix if present

for section in sections
    m = match(r"W(\d+)X(\d+)", section)
    if !isnothing(m)
        push!(depths, parse(Float64, m[1]))
        push!(weights, parse(Float64, m[2])) 
    end
end

closest_sections = String[]
for i in 1:lastindex(depths)
    closest_section = find_closest_section(depths[i], weights[i], available_sections)
    if isnothing(closest_section)
        println(row.depth, " ", row.weight)
    end
    push!(closest_sections, closest_section.name)
end

plot_slab(slab_params, closest_sections, text=true, mini=false, background=false, collinear=false, section_names=sections)