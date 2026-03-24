"""
    analyze_slab(self::SlabAnalysisParams; area::Real=-1., precision::Int64=4, axis::Union{Axis,Nothing}=nothing)

Analyzes the slab using the given cross sections based on the slab type. This function initializes a slab model and performs analysis based on the specified slab type, which can be isotropic, uniaxial, orthogonal, or orth_biaxial.

# Arguments
- `self::SlabAnalysisParams`: The parameters for slab analysis.
- `area::Real=-1.`: The area of the slab. If not specified, defaults to -1.
- `precision::Int64=4`: The precision for rounding results. Defaults to 4.

# Returns
- `self`: The updated `SlabAnalysisParams` object after analysis.

# Notes
- The function handles different slab types and calls the appropriate analysis function.
- It also manages NaN values in loads and finalizes the analysis by scaling if necessary.
"""
function analyze_slab(self::SlabAnalysisParams; area::Real=-1.)
    
    self.model.loads = Asap.AbstractLoad[]
    self.perp_vector_1d = [self.vector_1d[2], -self.vector_1d[1]]

    if self.plot_context.plot && self.plot_context.ax == nothing
        self.plot_context = setup_plot(self.plot_context.plot; size=(1200, 800))
    end

    if self.record_tributaries
        for element in self.model.elements[:beam]
            element_id = (element.nodeStart.nodeID, element.nodeEnd.nodeID)
            self.trib_dictionary["cycles"] = []
            self.trib_dictionary["interior_polygons"] = []
            self.trib_dictionary[element_id] = Dict(
                "trib_lines" => [], 
                "trib_area" => [], 
            )
        end
    end

    if self.slab_type == :isotropic
        self = analyze_isotropic_slab(self)
    elseif self.slab_type == :uniaxial
        self = analyze_uniaxial_slab(self)
    elseif self.slab_type == :orth_biaxial
        self = analyze_orth_biaxial_slab(self)    
    end

    if self.plot_context.plot
        plot_model(self.model, plot_context=self.plot_context)
        display(self.plot_context.fig)
    end

    handle_nan_loads(self.model.loads)
    Asap.solve!(self.model, reprocess=true)

    self = scale_analyze_slab(self, area)
    return self
end

function analyze_isotropic_slab(self)
    self = get_slab_loads_isotropic(self)
    self.load_dictionary = get_load_dictionary_by_id(self)
    self = consolidate_beam_loads(self)
    return self
end

function analyze_uniaxial_slab(self)
    self.perp = false
    self = get_slab_loads_uniaxial(self)
    self.load_dictionary = get_load_dictionary_by_id(self)
    self = consolidate_beam_loads(self)
    return self
end

function analyze_orth_biaxial_slab(self)
    println("Parallel vector: $(self.vector_1d)")
    self.perp = false
    self = get_slab_loads_uniaxial(self)
    display(self.plot_context.fig)
    println("Perpendicular vector: $(self.perp_vector_1d)")
    self.perp = true
    self = get_slab_loads_uniaxial(self)

    for load in self.model.loads
        load.value *= 0.5
    end

    self.load_areas .*= 0.5
    self.load_volumes .*= 0.5

    self.load_dictionary = get_load_dictionary_by_id(self)
    self = consolidate_beam_loads(self)
    
    return self
end

function consolidate_beam_loads(self::SlabAnalysisParams)
    tol = 1e-2
    loads = Asap.AbstractLoad[]
    elements = self.model.elements[:beam]

    for element in elements
        element_loads = get(self.load_dictionary, get_element_id(element), Asap.AbstractLoad[])
        if isempty(element_loads)
            continue
        end

        params, values = consolidate_loads(element_loads, tol)
        new_element_loads = [PointLoad(element, params[n], [0.0, 0.0, -values[n]]) for n in 1:lastindex(params)]
        self.load_dictionary[get_element_id(element)] = new_element_loads
    end

    all_loads = Asap.AbstractLoad[]
    for v in values(self.load_dictionary)
        append!(all_loads, v)
    end
    self.model.loads = all_loads

    return self
end

function consolidate_loads(beam_loads, tol)
    params = [beam_loads[1].position]
    values = [abs(beam_loads[1].value[3] / 2)]

    for j in 2:lastindex(beam_loads)
        close = false
        for k in 1:lastindex(params)
            if abs(beam_loads[j].position - params[k]) < tol
                close = true
                values[k] += abs(beam_loads[j].value[3] / 2)
                break
            end
        end
        if !close
            push!(params, beam_loads[j].position)
            push!(values, abs(beam_loads[j].value[3] / 2))
        end
    end

    return params, values
end

function handle_nan_loads(loads)
    for load in loads
        if isnan(load.value[3])
            load.value[3] = 0.0
        end
    end
end

function scale_analyze_slab(self, area)
    if area <= 0
        return self
    else
        println("Scaling required from $(self.area) m² to $area m².")
        if abs(area - self.area) > 1
            println("... scaling from $(self.area) m² to $area m² ...")
            self = scale_slab(self, area / self.area)
            self = analyze_slab(self)
        end
        return self
    end
end

function get_load_dictionary_by_element(self)
    load_dict = Dict{Element, Vector{Asap.AbstractLoad}}()

    for load in self.model.loads
        push!(get!(load_dict, load.element, Asap.AbstractLoad[]), load)
    end

    return load_dict
end

function get_element_id(element::Element)
    # Ensure the smaller node ID is first
    start_id = element.nodeStart.nodeID
    end_id = element.nodeEnd.nodeID
    return start_id < end_id ? (start_id, end_id) : (end_id, start_id)
end

function get_element_id(model::Asap.Model)
    element_id_lookup = Dict{Tuple{Int,Int}, Element}()
    for element in model.elements[:beam]
        element_id_lookup[get_element_id(element)] = element
    end
    return element_id_lookup
end

function get_element_from_id(self::SlabAnalysisParams, element_id::Tuple{Int,Int})
    return get(self.element_id_lookup, element_id, nothing)
end

function get_load_dictionary_by_id(self::SlabAnalysisParams)
    load_dict = Dict{Tuple{Int,Int}, Vector{Asap.AbstractLoad}}()

    for element in self.model.elements
        get!(load_dict, get_element_id(element), Asap.AbstractLoad[]) 
    end

    for load in self.model.loads
        element_id = get_element_id(load.element)
        push!(get!(load_dict, element_id, Asap.AbstractLoad[]), load)
    end

    return load_dict
end

function get_load_dictionary_by_id(model::Asap.Model)
    load_dict = Dict{Tuple{Int,Int}, Vector{Asap.AbstractLoad}}()

    for element in model.elements
        get!(load_dict, get_element_id(element), Asap.AbstractLoad[]) 
    end

    for load in model.loads
        element_id = get_element_id(load.element)
        push!(get!(load_dict, element_id, Asap.AbstractLoad[]), load)
    end

    return load_dict
end

function get_load_dictionary_by_index(self)
    load_dict = Dict{Int, Vector{Asap.AbstractLoad}}()

    element_to_index = Dict{Element, Int}()
    for (i, el) in enumerate(self.model.elements)
        element_to_index[el] = i
    end

    for load in self.model.loads
        idx = get(element_to_index, load.element, nothing)
        if idx !== nothing
            push!(get!(load_dict, idx, Asap.AbstractLoad[]), load)
        end
    end

    return load_dict
end
