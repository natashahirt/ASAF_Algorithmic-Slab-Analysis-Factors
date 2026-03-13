# Include necessary modules
include("../../SlabDesignFactors/scripts/_scripts.jl")

# Activate CairoMakie for plotting
CairoMakie.activate!()

# Define the path to the JSON file containing slab geometry
#path = "Geometries/special/topopt_test1.json"  # Update this path as needed
#path = "Geometries/nova/e4c4.json"

"""name = basename(splitext(path)[1])    # Name for the plot
#name = "topopt_test1_demo"
# Parse geometry from JSON
geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))
geometry, type_information = generate_from_json(geometry_dict, plot=true, drawn=false);"""

tiny = 1e-6
tiny_section = W_imperial("W0X0", tiny, tiny, tiny, tiny, tiny, tiny, tiny, tiny, tiny, tiny, tiny, tiny, tiny, tiny, tiny)
catalog_with_zero = vcat(allW_imperial(), [tiny_section])

begin

    h = 8
    w = 12
    resolution = 15
    geometry = generate_ground_structure(h, w, resolution, no_diagonal_connections=true, plot=true);
    name = "$(h)_$(w)_$(resolution)"

    slab_params = SlabAnalysisParams(
        geometry, 
        slab_name=name,
        slab_type=:isotropic,
        load_type=:indeterminate,
        vector_1d=[0,1], 
        slab_sizer=:uniform,
        spacing=.1, 
        plot_analysis=false,
        fix_param=true, 
        slab_units=:m,
    );

    # Sizing parameters
    beam_sizing_params = SlabSizingParams(
        live_load=psf_to_ksi(60), # ksi
        superimposed_dead_load=psf_to_ksi(20), # ksi
        live_factor=1.6, # -
        dead_factor=1.2, # -
        beam_sizer=:continuous,
        max_depth=40, # in
        beam_units=:in, # in, etc.
        serviceability_lim=360,
        collinear=false,
        minimum_continuous=false,
        n_max_sections=0,
        deflection_limit=false,
        catalog_discrete=catalog_with_zero
    );    

    # Create directory if it doesn't exist
    dir = "Documentation/Algorithms/Final_Algorithm_optimization_$(slab_params.slab_name)_optimizer"
    mkpath(dir)

    results_list = []
    slab_params_topopt = analyze_slab(slab_params)
    beam_sizing_params_topopt = beam_sizing_params

    starting_sections = ["W14X22" for _ in 1:length(slab_params.model.elements[:beam])]
    initial_vars = [get_geometry_vars(W_imperial(starting_sections[i])) for i in 1:length(starting_sections)]
    results_list = [initial_vars]
    starting_volume = get_volume(initial_vars, slab_params.model)
    starting_norm_mass = get_norm_mass(starting_volume, slab_params.area)
    fig = plot_slab(slab_params, starting_sections, text=false, mini=false, background=false, collinear=false)
    save("$dir/$(starting_norm_mass)_Starting_Section.pdf", fig)

    blend_ratio = 0.5

    for loop in 1:20

        for i in 1:20
            
            slab_params_topopt = calculate_slab_loads_indeterminate(slab_params_topopt)

            if length(results_list) == 1
                initial_vars = results_list[end]
            else
                initial_vars = [
                    blend_ratio .* results_list[end][j] .+ (1 - blend_ratio) .* results_list[end-1][j]
                    for j in 1:length(results_list[end])
                ]
            end

            slab_params_topopt, beam_sizing_params_topopt = optimal_beamsizer(slab_params_topopt, beam_sizing_params_topopt, initial_vars = initial_vars) # differentiable, uses Ipopt
            results = beam_sizing_params_topopt.minimizers

            # Apply section assignment once
            for (i, element) in enumerate(slab_params_topopt.model.elements[:beam])
                I_symm_section = I_symm(results[i]...)
                asap_section = Section(I_symm_section.A, steel_ksi.E, steel_ksi.G, I_symm_section.Ix, I_symm_section.Iy, I_symm_section.J)
                element.section = asap_section
            end

            push!(results_list, results)
            result_volume = get_volume(results, slab_params_topopt.model)
            result_norm_mass = get_norm_mass(result_volume, slab_params_topopt.area)
            println("Volume beams: $(result_volume)")
            println("Norm mass: $(result_norm_mass)")

            fig = plot_slab(slab_params_topopt, results_list[end], text=false, mini=false, background=false, collinear=false)
            save("$dir/$(result_norm_mass)_Iteration_$(i).pdf", fig)

            if i > 1    
                if abs(result_norm_mass - get_norm_mass(get_volume(results_list[i-1], slab_params_topopt.model), slab_params_topopt.area)) < .001
                    println("Converged at inner iteration $i")
                    break
                end
            end

            areas = beam_sizing_params_topopt.minimums
            println("Minimum area: $(minimum(areas))")
        
            println("--------------------")
        end

        volume_files = readdir(dir)

        minimizers, beam_sizing_params_topopt = optimize_indeterminate(slab_params_topopt, beam_sizing_params_topopt, initial_vars = results_list[end]);
        
        push!(results_list, minimizers)
        
        result_volume = get_volume(results_list[end], slab_params_topopt.model)
        result_norm_mass = get_norm_mass(result_volume, slab_params_topopt.area)
        println("Volume beams: $(result_volume)")
        println("Norm mass: $(result_norm_mass)")
        areas = [I_symm(minimizers[i]...).A for i in 1:lastindex(minimizers)]
        println("Minimum area: $(minimum(areas))")
        fig = plot_slab(slab_params_topopt, minimizers, text=false, mini=false, background=false, collinear=false)
        save("$dir/$(result_norm_mass)_Topopt_$(loop).pdf", fig)

        # Check for convergence with best result so far in directory        
        if !isempty(volume_files)
            # Filter out .DS_Store and other hidden files, only keep files with numeric volume prefix
            volumes = map(f -> parse(Float64, split(f, "_")[1]), filter(f -> !startswith(f, "."), volume_files))
            best_volume = minimum(volumes)
            best_norm_mass = get_norm_mass(best_volume, slab_params_topopt.area)
            
            # Compare current volume with best
            if abs(result_norm_mass - best_norm_mass) < 0.001
                println("Converged with best result ($(best_norm_mass)) after $loop outer loop iterations.")
                break
            end
        end

    end

end

h = 8
w = 12
resolution = 15
geometry = generate_ground_structure(h, w, resolution, no_diagonal_connections=false, plot=true, moment_release=true);

begin

    h = 8
    w = 12
    resolution = 15
    geometry = generate_ground_structure(h, w, resolution, no_diagonal_connections=false, plot=true, moment_release=true);
    name = "$(h)_$(w)_$(resolution)"

    slab_params = SlabAnalysisParams(
        geometry, 
        slab_name=name,
        slab_type=:isotropic,
        load_type=:indeterminate,
        vector_1d=[0,1], 
        slab_sizer=:uniform,
        spacing=.1, 
        plot_analysis=false,
        fix_param=true, 
        slab_units=:m,
    );

    # Sizing parameters
    beam_sizing_params = SlabSizingParams(
        live_load=psf_to_ksi(60), # ksi
        superimposed_dead_load=psf_to_ksi(20), # ksi
        live_factor=1.6, # -
        dead_factor=1.2, # -
        beam_sizer=:continuous,
        max_depth=40, # in
        beam_units=:in, # in, etc.
        serviceability_lim=360,
        collinear=false,
        minimum_continuous=false,
        n_max_sections=0,
        deflection_limit=false,
        catalog_discrete=catalog_with_zero
    );    

    # Create directory if it doesn't exist
    dir = "Documentation/Algorithms/Final_Algorithm_optimization_$(slab_params.slab_name)_drawn_diag_asymm_moment_release"
    mkpath(dir)

    results_list = []
    slab_params_topopt = analyze_slab(slab_params)
    beam_sizing_params_topopt = beam_sizing_params

    starting_sections = ["W14X22" for _ in 1:length(slab_params.model.elements[:beam])]
    initial_vars = [get_geometry_vars(W_imperial(starting_sections[i])) for i in 1:length(starting_sections)]
    results_list = [initial_vars]
    starting_volume = get_volume(initial_vars, slab_params.model)
    starting_norm_mass = get_norm_mass(starting_volume, slab_params.area)
    fig = plot_slab(slab_params, starting_sections, text=false, mini=false, background=false, collinear=false)
    save("$dir/$(starting_norm_mass)_Starting_Section.pdf", fig)

    blend_ratio = 0.5

    for loop in 1:20

        for i in 1:20
            
            slab_params_topopt = calculate_slab_loads_indeterminate(slab_params_topopt)

            if length(results_list) == 1
                initial_vars = results_list[end]
            else
                initial_vars = [
                    blend_ratio .* results_list[end][j] .+ (1 - blend_ratio) .* results_list[end-1][j]
                    for j in 1:length(results_list[end])
                ]
            end

            """if i % 2 == 0
                beam_sizing_params.beam_sizer = :continuous
            else
                beam_sizing_params.beam_sizer = :discrete
            end"""

            slab_params_topopt, beam_sizing_params_topopt = optimal_beamsizer(slab_params_topopt, beam_sizing_params_topopt, initial_vars = initial_vars) # differentiable, uses Ipopt
            results = beam_sizing_params_topopt.minimizers

            # Apply section assignment once
            for (i, element) in enumerate(slab_params_topopt.model.elements[:beam])
                I_symm_section = I_symm(results[i]...)
                asap_section = Section(I_symm_section.A, steel_ksi.E, steel_ksi.G, I_symm_section.Ix, I_symm_section.Iy, I_symm_section.J)
                element.section = asap_section
            end

            push!(results_list, results)
            result_volume = get_volume(results, slab_params_topopt.model)
            result_norm_mass = get_norm_mass(result_volume, slab_params_topopt.area)
            println("Volume beams: $(result_volume)")
            println("Norm mass: $(result_norm_mass)")

            fig = plot_slab(slab_params_topopt, results_list[end], text=false, mini=false, background=false, collinear=false)
            save("$dir/$(result_norm_mass)_Iteration_$(i).pdf", fig)

            if i > 1    
                if abs(result_norm_mass - get_norm_mass(get_volume(results_list[i-1], slab_params_topopt.model), slab_params_topopt.area)) < .001
                    println("Converged at inner iteration $i")
                    break
                end
            end

            areas = beam_sizing_params_topopt.minimums
            println("Minimum area: $(minimum(areas))")
        
            println("--------------------")
        end

        copy_results = deepcopy(results_list[end])

        # Randomly perturb each variable in minimizers by a different amount
        for i in 1:length(copy_results)
            copy_results[i][1] += (1*rand()) # Randomly increase h by up to 10%
            copy_results[i][2] += (1*rand()) # Randomly increase w by up to 10% 
            copy_results[i][3] += (1*rand()) # Randomly increase tw by up to 10%
            copy_results[i][4] += (1*rand()) # Randomly increase tf by up to 10%
        end

        push!(results_list, copy_results)

    end

end


# Define paths to result directories
coarse_grid_path = "Documentation/Algorithms/Final_Algorithm_optimization_8_12_5_grid"
coarse_diag_path = "Documentation/Algorithms/Final_Algorithm_optimization_8_12_5_convergence_optim"
fine_grid_path = "Documentation/Algorithms/Final_Algorithm_optimization_8_12_15_optimizer"
fine_diag_path = "Documentation/Algorithms/Final_Algorithm_optimization_8_12_15_random_shuffle_diag"

# Function to extract volume from filename and sort by creation time
function get_sorted_volumes(path; normalized::Bool=true)
    files = readdir(path)
    files = filter(f -> !startswith(f, "."), files)
    
    # Get creation times and volumes
    file_data = map(files) do f
        ctime = stat(joinpath(path, f)).ctime
        vol = parse(Float64, split(f, "_")[1])
        (ctime, vol)
    end
    
    # Sort by creation time and extract volumes
    sort!(file_data, by=x->x[1])
    volumes = [x[2] for x in file_data]
    
    # Normalize to percentages
    if normalized == true
        volumes = volumes ./ volumes[1] .* 100
    end
    
    # Return volumes up to convergence
    return volumes[1:21]
end

# Create vectors of volumes sorted by creation time
coarse_grid_volumes = get_sorted_volumes(coarse_grid_path, normalized=true)
coarse_diag_volumes = get_sorted_volumes(coarse_diag_path, normalized=true) 
fine_grid_volumes = get_sorted_volumes(fine_grid_path, normalized=true)
fine_diag_volumes = get_sorted_volumes(fine_diag_path, normalized=true)

max_volumes = [maximum(coarse_grid_volumes), maximum(coarse_diag_volumes), maximum(fine_grid_volumes), maximum(fine_diag_volumes)]
println(max_volumes)

fig = Figure(size=(1000, 600));
ax = Axis(fig[1, 1], xlabel="Iteration #", ylabel="% Decrease in structural mass");

scatterlines!(ax, 0:20, coarse_grid_volumes, label="Coarse Grid", color=色[:ceruleanblue])
scatterlines!(ax, 0:20, coarse_diag_volumes, label="Coarse Diagonal", color=色[:irispurple]) 
scatterlines!(ax, 0:20, fine_grid_volumes, label="Fine Grid", color=色[:powderblue])
scatterlines!(ax, 0:20, fine_diag_volumes, label="Fine Diagonal", color=色[:magenta])

axislegend(ax, framecolor=:white)
xlims!(ax, 0, 20)

fig
save("SlabDesignFactors/plot/figures/convergence_plots.pdf", fig)


# Get unique coordinates from raster_df centerpoint_coords
slab_params = analyze_slab(slab_params);
get_raster_df(slab_params; resolution=50)
"""process_continuous_beams_topopt(slab_params, beam_sizing_params)

design_variables = vcat([get_geometry_vars(W_imperial("W8X35")) for _ in 1:length(slab_params.model.elements[:beam])]...)

for (i,element) in enumerate(slab_params.model.elements[:beam])
    I_symm_section = I_symm(design_variables[4*i-3:4*i]...)
    asap_section = Section(I_symm_section.A, steel_ksi.E, steel_ksi.G, I_symm_section.Ix, I_symm_section.Iy, I_symm_section.J)
    element.section = asap_section
end

slab_params, beam_sizing_params = optimal_beamsizer(slab_params, beam_sizing_params);
results = postprocess_slab(slab_params, beam_sizing_params);
objective_variables = [results.embodied_carbon_slab, results.embodied_carbon_beams]"""

# REDISTRIBUTION

slab_params = analyze_slab(slab_params);
slab_params, beam_sizing_params = optimal_beamsizer(slab_params, beam_sizing_params);
results_standard_optimization = postprocess_slab(slab_params, beam_sizing_params);
results_standard_optimization.embodied_carbon_beam
centerpoint_coords = identify_raster_points(slab_params, 5)
scatter(Point2f.(centerpoint_coords), color=:red, markersize=8, axis=(limits=(0, 12, 0, 10),))

slab_params_topopt, beam_sizing_params_topopt = optimize_indeterminate(slab_params, beam_sizing_params, initial_vars = []);
fig = plot_slab(slab_params_topopt, beam_sizing_params_topopt, text=false, mini=false, background=false, collinear=false)

for result in results_list
    volume = get_volume(result.minimizers, slab_params.model)
    println(volume)
end


function optimize_beam_sections(slab_params, beam_sizing_params)

    slab_params = analyze_slab(slab_params);

    model = slab_params.model
    beam_elements = model.elements[:beam]
    default_section = W_imperial("W8X35")
    
    design_variables = vcat([get_geometry_vars(default_section) for _ in 1:length(beam_elements)]...)

    function objective(design_variables)
        i_sections = [I_symm(design_variables[4*i-3:4*i]...) for i in 1:length(beam_elements)]
        
        for (i, element) in enumerate(slab_params.model.elements[:beam])
            section = Section(i_sections[i].A, steel_ksi.E, steel_ksi.G, i_sections[i].Ix, i_sections[i].Iy, i_sections[i].J)
            
            # Ensure non-negative area and moments of inertia
            if section.A < 0 || section.Ix < 0 || section.Iy < 0 || section.J < 0
                println("Warning: Negative area or moment of inertia detected. Applying penalty.")
                return Inf  # Return a large penalty value
            end
            
            element.section = section
        end

        slab_params = calculate_slab_loads_indeterminate(slab_params)
        slab_params, beam_sizing_params = optimal_beamsizer(slab_params, beam_sizing_params, initial_vars = [design_variables[4*i-3:4*i] for i in 1:length(beam_elements)])
        
        A = [A_I_symm(vars...) for vars in beam_sizing_params.minimizers]
        
        total_volume = sum(A[i] * beam_elements[i].length for i in 1:length(beam_elements))
        
        return total_volume
    end

    if beam_sizing_params.minimum_continuous == true
        min_h, min_w, min_tw, min_tf = get_geometry_vars(W_imperial("W6X8.5"))
    else
        min_h, min_w, min_tw, min_tf = [0.01, 0.01, 0.001, 0.001]
    end

    max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))

    bounds = [(min_h, max_h), (min_w, max_w), (min_tw, max_tw), (min_tf, max_tf)]

    search_range = [(Float64(min), Float64(max)) for (min, max) in repeat(bounds, length(beam_elements))]

    result = bboptimize(objective, SearchRange=search_range, MaxSteps=3600, Method=:adaptive_de_rand_1_bin_radiuslimited)

    return beam_sizing_params

end

using BlackBoxOptim
beam_sizing_params_topopt = optimize_beam_sections(slab_params, beam_sizing_params);

results_topopt = postprocess_slab(slab_params_topopt, beam_sizing_params_topopt);
results_topopt.embodied_carbon_beams

results_topopt.minimizers
print_forces(results_topopt)

areas = sum([I_symm(results_topopt.minimizers[i]...).A for i in 1:n_beams])

fig = plot_slab(slab_params_topopt, beam_sizing_params_topopt, text=false, mini=false, background=false, collinear=false)

volumes = sum([I_symm(beam_sizing_params_topopt.minimizers[i]...).A * beam_sizing_params_topopt.model.elements[:beam][i].length * 1/convert_to_m[:in] for i in 1:n_beams])

sum(results_areas)
sum(results_parallel_areas)
sum(results_2_areas)

slab_params = analyze_slab(slab_params);
slab_params, beam_sizing_params = optimal_beamsizer(slab_params, beam_sizing_params, initial_vars = results_topopt.minimizers);
beam_sizing_params_volumes = [I_symm(beam_sizing_params.minimizers[i]...).A * beam_sizing_params.model.elements[:beam][i].length for i in 1:n_beams]
sum(beam_sizing_params_volumes)

results = postprocess_slab(slab_params, beam_sizing_params);
volume_optimized = sum([I_symm(results.minimizers[i]...).A * beam_sizing_params.model.elements[:beam][i].length * 1/convert_to_m[:in] for i in 1:n_beams])
sum([I_symm(results.minimizers[i]...).A for i in 1:n_beams])
print_forces(results)
results.minimizers
results.embodied_carbon_beams

fig = plot_slab(slab_params, beam_sizing_params, text=false, mini=false, background=false, collinear=false)
