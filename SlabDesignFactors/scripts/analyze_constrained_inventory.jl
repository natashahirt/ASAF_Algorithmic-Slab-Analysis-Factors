# Include necessary modules
include("_scripts.jl")

# Activate CairoMakie for plotting
CairoMakie.activate!()

live_load_dict = Dict(
        "office" => 60,
        "school" => 50,
        "warehouse" => 250,
    )

# Define the path to the JSON file containing slab geometry
# path = "Geometries/special/triple_bay.json"  # Update this path as needed

for name in ["office", "school", "warehouse"]

    path = "Geometries/validation/$(name).json"
    name = basename(splitext(path)[1])    # Name for the plot
    # Parse geometry from JSON
    geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict));
    geometry, type_information = generate_from_json(geometry_dict, plot=true, drawn=true);

    live_load_dict = Dict(
        "office" => 60,
        "school" => 50,
        "warehouse" => 250,
    )

    for i in range(0, 30)

        # Analyze the slab to get dimensions
        slab_params = SlabAnalysisParams(
            geometry, 
            slab_name=name,
            slab_type=:uniaxial,
            slab_thickness=4.75 * convert_to_m[:in], # Convert 4.75 inches to meters
            vector_1d=[0,1], 
            slab_sizer=:uniform,
            spacing=.1, 
            plot_analysis=true,
            fix_param=true, 
            slab_units=:m,
            i_holes=type_information["i_holes"],
            i_perimeter=type_information["i_perimeter"],
        );

        # Sizing parameters
        beam_sizing_params = SlabSizingParams(
            live_load=psf_to_ksi(live_load_dict[name]),
            superimposed_dead_load=psf_to_ksi(20), # ksi
            slab_dead_load=psf_to_ksi(45), # ksi
            façade_load=plf_to_kpi(500), # kip/in
            live_factor=1.6, # -
            dead_factor=1.2, # -
            beam_sizer=:discrete,
            max_depth=40, # in
            beam_units=:in, # in, etc.
            serviceability_lim=360,
            collinear=false,
            drawn=true,
            element_ids=type_information["element_ids"],
            minimum_continuous=true,
            n_max_sections=0,
        );

        slab_params = analyze_slab(slab_params);

        beam_sizing_params.n_max_sections = i

        slab_params_analyzed, beam_sizing_params_analyzed = optimal_beamsizer(slab_params, beam_sizing_params);
        
        if isempty(beam_sizing_params_analyzed.minimizers)
            continue
        end

        slab_results = postprocess_slab(slab_params_analyzed, beam_sizing_params_analyzed);
        mkpath("SlabDesignFactors/results/constrained_inventory_resized/")
        append_results_to_csv("SlabDesignFactors/results/constrained_inventory_resized/", "$name", [slab_results], unique_sections=i)

        fig = plot_slab(slab_params_analyzed, beam_sizing_params_analyzed, text=true, mini=false, background=false, collinear=false)
        display(fig)
        println(slab_results.norm_mass_beams)

        # Create directory if it doesn't exist
        mkpath("Documentation/constrained_inventory_resized/$(name)")
        save("Documentation/constrained_inventory_resized/$(name)/$(name)_$(i).pdf", fig)
        
        unique_sections = unique(beam_sizing_params_analyzed.minimizers)
        println("Number of unique sections: $(length(unique_sections))")

    end

end


# Maximum depth
for name in ["office", "school", "warehouse"]

    path = "Geometries/validation/$(name).json"
    name = basename(splitext(path)[1])    # Name for the plot
    # Parse geometry from JSON
    geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict));
    geometry, type_information = generate_from_json(geometry_dict, plot=true, drawn=true);

    live_load_dict = Dict(
        "school" => 50,
        "office" => 60,
        "warehouse" => 250,
    )

    for i in [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 1000]

        # Analyze the slab to get dimensions
        slab_params = SlabAnalysisParams(
            geometry, 
            slab_name=name,
            slab_type=:uniaxial,
            slab_thickness=4.75 * convert_to_m[:in], # Convert 4.75 inches to meters
            vector_1d=[0,1], 
            slab_sizer=:uniform,
            spacing=.1, 
            plot_analysis=true,
            fix_param=true, 
            slab_units=:m,
            i_holes=type_information["i_holes"],
            i_perimeter=type_information["i_perimeter"],
        );

        # Sizing parameters
        beam_sizing_params = SlabSizingParams(
            live_load=psf_to_ksi(live_load_dict[name]),
            superimposed_dead_load=psf_to_ksi(20), # ksi
            slab_dead_load=psf_to_ksi(45), # ksi
            façade_load=plf_to_kpi(500), # kip/in
            live_factor=1.6, # -
            dead_factor=1.2, # -
            beam_sizer=:discrete,
            max_depth=i, # in
            beam_units=:in, # in, etc.
            serviceability_lim=360,
            collinear=false,
            drawn=true,
            element_ids=type_information["element_ids"],
            minimum_continuous=true,
            n_max_sections=0,
        );

        slab_params = analyze_slab(slab_params);
        slab_params, beam_sizing_params = optimal_beamsizer(slab_params, beam_sizing_params);
        
        if isempty(beam_sizing_params.minimizers)
            continue
        end

        slab_results = postprocess_slab(slab_params, beam_sizing_params);
        mkpath("SlabDesignFactors/results/constrained_inventory_max_depth/")
        append_results_to_csv("SlabDesignFactors/results/constrained_inventory_max_depth/", "$name", [slab_results], unique_sections=i)

        fig = plot_slab(slab_params, beam_sizing_params, text=true, mini=false, background=false, collinear=false)
        # Create directory if it doesn't exist
        mkpath("Documentation/constrained_inventory_max_depth/$(name)")
        save("Documentation/constrained_inventory_max_depth/$(name)/$(name)_$(i).pdf", fig)
        
        unique_sections = unique(beam_sizing_params.minimizers)
        println("Number of unique sections: $(length(unique_sections))")

    end

end


name = "warehouse"
path = "Geometries/validation/$(name).json"
# Parse geometry from JSON
geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict));
geometry, type_information = generate_from_json(geometry_dict, plot=true, drawn=true);

slab_params = SlabAnalysisParams(
    geometry, 
    slab_name=name,
    slab_type=:uniaxial,
    vector_1d=[0,1], 
    slab_thickness=4.75 * convert_to_m[:in], # Convert 4.75 inches to meters
    slab_sizer=:uniform,
    spacing=.1, 
    plot_analysis=true,
    fix_param=true, 
    slab_units=:m,
    i_holes=type_information["i_holes"],
    i_perimeter=type_information["i_perimeter"],
);

# Sizing parameters
beam_sizing_params = SlabSizingParams(
    live_load=psf_to_ksi(live_load_dict[name]),
    superimposed_dead_load=psf_to_ksi(20), # ksi
    slab_dead_load=psf_to_ksi(45), # ksi
    façade_load=plf_to_kpi(500), # kip/in
    live_factor=1.6, # -
    dead_factor=1.2, # -
    beam_sizer=:discrete,
    max_depth=1000, # in
    beam_units=:in, # in, etc.
    serviceability_lim=360,
    collinear=false,
    drawn=true,
    element_ids=type_information["element_ids"],
    minimum_continuous=true,
    n_max_sections=0,
);

slab_params = analyze_slab(slab_params);
slab_params, beam_sizing_params = optimal_beamsizer(slab_params, beam_sizing_params);
slab_results = postprocess_slab(slab_params, beam_sizing_params);
print_forces(slab_results)

fig = plot_slab(slab_params, beam_sizing_params, text=true, mini=false, background=false, collinear=false)
#save("SlabDesignFactors/plot/figures/constrained_inventory/office_numbered.pdf", fig)

#fig = plot_slab(slab_params, beam_sizing_params, text=true, mini=false, background=false, collinear=false)
#save("SlabDesignFactors/plot/figures/constrained_inventory/school_numbered.pdf", fig)

source_path = "Geometries/validation/$(name)_original.csv"
df = CSV.read(source_path, DataFrame, header=false)
rename!(df, [:depth, :weight])

# Create DataFrame with depth and weight columns
df = DataFrame(depth=df.depth, weight=df.weight)

total_weight = sum(df.weight .* [element.length * 3.28084 for element in slab_params.model.elements[:beam]])
total_weight_kg = total_weight * 0.453592
total_weight_norm = total_weight_kg / 2297.9

# Convert depth and weight columns to W-section format
df.section = map(row -> "W$(row.depth)X$(trunc(Int, row.weight))", eachrow(df))
# Get all available W sections
available_sections = allW_imperial()

closest_sections = String[]
# Function to find closest section by depth and weight
function find_closest_section(target_depth, target_weight, available_sections)
    min_diff = Inf
    closest_section = nothing
    
    for section in available_sections
        # Extract depth and weight from section name
        m = match(r"W(\d+)X(\d+)", string(section.name))
        if isnothing(m)
            continue
        end
        depth = parse(Float64, m[1])
        weight = parse(Float64, m[2])
        
        # Calculate difference using depth and weight
        diff = sqrt((depth - target_depth)^2 + (weight - target_weight)^2)
        
        if diff < min_diff
            min_diff = diff
            closest_section = section
        end
    end
    return closest_section
end

closest_sections = String[]
for row in eachrow(df)
    closest_section = find_closest_section(row.depth, row.weight, available_sections)
    if isnothing(closest_section)
        println(row.depth, " ", row.weight)
    end
    push!(closest_sections, closest_section.name)
end

fig = plot_slab(slab_params, closest_sections, text=true, mini=false, background=false, collinear=false, section_names=df.section)
save("SlabDesignFactors/plot/figures/constrained_inventory/$(name)_original.pdf", fig)

area_analysis = [W_imperial(section).A * convert_to_m[:in]^2 for section in slab_results.sections]
area_original = [W_imperial(section).A * convert_to_m[:in]^2 for section in closest_sections]

volume_analysis = area_analysis .* [element.length for element in slab_params.model.elements[:beam]]
volume_original = area_original .* [element.length for element in slab_params.model.elements[:beam]]

mass_analysis = volume_analysis .* ρ_STEEL
mass_original = volume_original .* ρ_STEEL

norm_mass_analysis = mass_analysis ./ slab_params.area
norm_mass_original = mass_original ./ slab_params.area

# Calculate percentage difference in mass
sum(norm_mass_analysis)
sum(norm_mass_original)
norm_mass_percent_diff = (sum(norm_mass_analysis) - sum(norm_mass_original)) / sum(norm_mass_original) * 100

norm_mass_delta = norm_mass_analysis .- norm_mass_original
area_delta = area_analysis .- area_original

# Find index of maximum area delta
max_delta_idx = argmax(abs.(norm_mass_delta))
norm_mass_delta[max_delta_idx]

# Get the original and analysis sections at that index
original_section = closest_sections[max_delta_idx]
analysis_section = slab_results.sections[max_delta_idx]

# Count number of unique sections in original sections
n_unique_original = length(unique(closest_sections))
n_unique_analysis = length(unique(slab_results.sections))

fig = plot_slab_section_delta(slab_params, norm_mass_delta, text=true, mini=false, background=false, collinear=false)
save("SlabDesignFactors/plot/figures/constrained_inventory/$(name)_norm_mass_delta.pdf", fig)