# analyze_single_slab.jl
include("../SlabDesignFactors.jl")

using Base.Threads
using JSON
using CSV
using DataFrames
using .SlabDesignFactors

function analyze_experiments(results_path::String, completion_file::String)

    println("Dependencies loaded successfully.")

    # Read parameters from file
    json_path = "/home/nhirt/2024_Slab-Design-Factors/SlabDesignFactors/jsons/topology/"

    experiment_max_depths = Dict(
        :results_name => "max_depths",
        :slab_types => [:isotropic],    # Slab types
        :vector_1ds => [[0.0, 0.0]],    # Vectors
        :max_depths => [10, 15, 20, 25, 30, 35, 40, 45, 50, 10000],
        :slab_sizers => [:uniform, :cellular],
        :fix_params => [true]
    )

    experiment_fix_params = Dict(
        :results_name => "fix_params",
        :slab_types => [:isotropic, :orth_biaxial, :orth_biaxial, :uniaxial, :uniaxial, :uniaxial, :uniaxial],        # Slab types
        :vector_1ds => [[0.0, 0.0], [1.0, 0.0,], [1.0, 1.0,], [1.0, 0.0,], [0.0, 1.0,], [1.0, 1.0,], [1.0, -1.0,]],    # Vectors
        :max_depths => [25, 40],
        :slab_sizers => [:uniform, :cellular],
        :fix_params => [false]
    )

    experiment_dicts = [experiment_max_depths]

    # Function to process a single set of parameters
    function process_params(json_path, results_path, experiment_dict)

        # Get the JSON path and results name
        results_name = experiment_dict[:results_name]
        slab_types = experiment_dict[:slab_types]
        vector_1ds = experiment_dict[:vector_1ds]
        max_depths = experiment_dict[:max_depths]
        slab_sizers = experiment_dict[:slab_sizers]
        fix_params = experiment_dict[:fix_params]

        println("Analyzing all JSONs in $(json_path)")
        println("Saving results to $(results_path)")
        println("Experiment name: $(results_name)")

        # Define the path to the JSON file containing slab geometry
        sub_paths = filter(x -> endswith(x, ".json"), readdir(json_path))
        println(sub_paths)
        
        # Evaluate slabs
        for max_depth in max_depths

            for slab_sizer in slab_sizers

                for fix_param in fix_params

                        for (i, slab_type) in enumerate(slab_types)

                        vector_1d = vector_1ds[i]

                        # Check if this configuration already exists in results file
                        results_file = results_path * results_name * ".csv"

                        for sub_path in sub_paths

                            path = json_path * sub_path
                            name = replace(sub_path, ".json" => "")  

                            if isfile(results_file)
                                existing_df = CSV.read(results_file, DataFrame)
                                if any((existing_df.name .== name) .&
                                        (existing_df.max_depth .== max_depth) .&
                                        (existing_df.slab_sizer .== string(slab_sizer)) .& 
                                        (existing_df.slab_type .== string(slab_type)) .&
                                        (existing_df.vector_1d_x .== vector_1d[1]) .&
                                        (existing_df.vector_1d_y .== vector_1d[2]))
                                    println("Already analyzed $(name) for $(slab_type) $(vector_1d) $(slab_sizer) $(max_depth) in.")
                                    continue
                                end
                            end

                            println("================================================")
                            println("$(name): $(slab_type) $(vector_1d) $(slab_sizer) $(max_depth)in")
                            println("================================================")

                            # Parse geometry from JSON
                            json_string = replace(read(path, String), "\\n" => "")
                            geometry_dict = JSON.parse(JSON.parse(json_string, dicttype=Dict))

                            # Use the function from the module
                            geometry = SlabDesignFactors.generate_from_json(geometry_dict; plot=false, drawn=false)

                            # Create and analyze the slab
                            slab_params = SlabDesignFactors.SlabAnalysisParams(
                                geometry, 
                                slab_name=name,
                                slab_type=slab_type,
                                vector_1d=vector_1d, 
                                slab_sizer=slab_sizer,
                                spacing=.1, 
                                plot_analysis=false,
                                fix_param=true, 
                                slab_units=:m,
                            );

                            # Sizing parameters
                            beam_sizing_params = SlabDesignFactors.SlabSizingParams(
                                live_load=SlabDesignFactors.psf_to_ksi(50), # ksi
                                superimposed_dead_load=SlabDesignFactors.psf_to_ksi(15), # ksi
                                live_factor=1.6, # -
                                dead_factor=1.2, # -
                                beam_sizer=:discrete, # iteration runs through both discrete and continuous
                                max_depth=max_depth, # in
                                beam_units=:in, # in, etc.
                                serviceability_lim=360,
                                minimum_continuous=true
                            );

                            iteration_result = collect(SlabDesignFactors.iterate_discrete_continuous(slab_params, beam_sizing_params));
                        
                            SlabDesignFactors.append_results_to_csv(results_path, String(results_name), iteration_result)

                            GC.gc() # garbage collect

                        end

                    end

                end

            end

        end

    end

    Threads.@threads for experiment_dict in experiment_dicts
        process_params(json_path, results_path, experiment_dict)
    end

    # Create a completion file to signal the end of processing
    open(completion_file, "w") do f
        write(f, "Experiments complete")
    end

end

# Main execution
function main()
    args = ARGS

    if length(args) != 2
        println("Usage: julia executable_experiments.jl <results_path> <completion_file>")
        return
    end

    results_path = args[1]
    completion_file = args[2]
    analyze_experiments(results_path, completion_file)

end

main()