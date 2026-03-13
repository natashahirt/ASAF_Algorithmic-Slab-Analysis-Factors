function get_tributary_lines(msg; tol=0.001)
    geometry_dict = JSON.parse(replace(msg, "\\n" => ""), dicttype=Dict)
    geometry = generate_from_json(geometry_dict, plot=false, drawn=false)
    name = "test"
    
    # Analyze the slab to get dimensions
    slab_params = SlabAnalysisParams(
        geometry, 
        slab_name=name,
        slab_type=:isotropic,
        vector_1d=[0,0], 
        slab_sizer=:uniform,
        spacing=geometry_dict["s"][1], 
        plot_analysis=true,
        fix_param=true, 
        slab_units=:m,
        record_tributaries=true,
    )

    slab_params = analyze_slab(slab_params)

    node_dictionary = Dict([(i, node.position) for (i, node) in enumerate(slab_params.model.nodes) if node.id != :fixed])

    geometry_dictionary = Dict(
        "nodes" => node_dictionary,
        "tributaries" => slab_params.trib_dictionary
    )

    return JSON.json(geometry_dictionary)
end