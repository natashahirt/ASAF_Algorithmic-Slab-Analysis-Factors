function get_max_deflection(test_result::Union{DataFrameRow, DataFrameRow{DataFrame, DataFrames.Index}}; plot_analysis::Bool=true)

    # Parse geometry from JSON
    path = "Geometries/$(test_result.category)/$(test_result.name).json"
    geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict));
    material = material_dict[:in]
    sections = parse_sections(String(test_result.sections))
    if !contains(sections[1], "W")
        println(sections[1])
        return [NaN], [NaN]
    end
    sections = [toASAPframe_W(W_imperial(W_dict_imperial[section]), material.E, material.G; convert=false) for section in sections]
    geometry = generate_from_json(geometry_dict, plot=false, drawn=false, sections=sections);

    # Analyze the slab to get dimensions
    slab_params = SlabAnalysisParams(
        geometry, 
        slab_name=test_result.name,
        slab_type=Symbol(test_result.slab_type),
        vector_1d=[test_result.vector_1d_x, test_result.vector_1d_y], 
        slab_sizer=Symbol(test_result.slab_sizer),
        spacing=.1, 
        plot_analysis=plot_analysis,
        fix_param=true, 
        slab_units=:m,
    );

    beam_sizing_params = SlabSizingParams(
        live_load=psf_to_ksi(50), # ksi
        superimposed_dead_load=psf_to_ksi(15), # ksi
        live_factor=1.6, # -
        dead_factor=1.2, # -
        beam_sizer=Symbol(test_result.beam_sizer),
        max_depth=Int64(test_result.max_depth), # in
        beam_units=:in, # in, etc.
        serviceability_lim=360,
        collinear=Bool(test_result.collinear),
        minimum_continuous=true
    );

    slab_params = analyze_slab(slab_params);
    
    # convert model lengths to inches
    conversion_factor = convert_to_m[slab_params.slab_units] * 1/convert_to_m[beam_sizing_params.beam_units]
    beam_sizing_params.area = slab_params.area * conversion_factor^2

    # Adjust max depth based on slab depth
    if beam_sizing_params.max_assembly_depth
        slab_depth = maximum(slab_params.slab_depths) * conversion_factor # Convert to inches
        beam_sizing_params.max_beam_depth = beam_sizing_params.max_depth - slab_depth
    end

    beam_sizing_params.model = get_scaled_model(slab_params, beam_sizing_params, conversion_factor);
    beam_sizing_params.load_dictionary = get_load_dictionary_by_id(beam_sizing_params.model);

    model = beam_sizing_params.model
    load_dictionary = beam_sizing_params.load_dictionary
    resolution = 200

    local_deflections = []
    global_deflections = []
    for element in model.elements[:beam]
        displacements = ElementDisplacements(element, load_dictionary[get_element_id(element)], resolution=resolution)
        ulocal = displacements.ulocal
        uglobal = displacements.uglobal
        push!(local_deflections, maximum(abs.(ulocal)))
        push!(global_deflections, maximum(abs.(uglobal)))
    end
    
    return local_deflections, global_deflections

end

function get_max_deflection(df::DataFrame)

    max_local_deflections = []
    max_global_deflections = []

    for row in eachrow(df)
        if row.area == 0 || isnan(row.area)
            push!(max_local_deflections, NaN)
            push!(max_global_deflections, NaN)
            continue
        end
        local_deflections, global_deflections = get_max_deflection(row, plot_analysis=false)
        if length(local_deflections) == 1
            continue
        end
        push!(max_local_deflections, maximum(local_deflections))
        push!(max_global_deflections, maximum(global_deflections))
    end

    return max_local_deflections, max_global_deflections

end