abstract type AbstractVariable end
abstract type IndependentVariable <: AbstractVariable end
abstract type AbstractOptParams end

"""
    PlotContext

Context for plotting.

# Fields
"""
mutable struct PlotContext
    plot::Bool                              # Plotting flag
    fig::Union{Figure,Nothing}              # Figure handle
    ax::Union{Axis,Nothing}                 # Axis handle

    function PlotContext(plot::Bool, fig::Union{Figure,Nothing}, ax::Union{Axis,Nothing})
        new(plot, fig, ax)
    end
end

"""
    SlabAnalysisParams

Parameters for slab analysis, including pre-sizing and post-sizing attributes.
"""
mutable struct SlabAnalysisParams <: AbstractOptParams
    model::Asap.Model               # The model
    slab_name::String               # Slab name
    slab_type::Symbol              # :uniaxial or :isotropic
    load_type::Symbol              # :determinate or :indeterminate
    vector_1d::Vector{<:Real}      # Direction of uniaxial slab
    slab_thickness::Real          # Slab thickness [in]
    perp::Bool                    # Whether to use the perpendicular direction
    perp_vector_1d::Vector{<:Real} # Perpendicular direction of orth_biaxial slab
    slab_sizer::Symbol             # :cellular or :uniform
    fix_param::Bool                # Adjust for stiffness of elements or stay at 0.5?
    spacing::Float64               # Density of sampling [-]
    area::Float64                  # Sum area of slab [-²]
    areas::Vector{<:Real}         # Areas of each cell [-²]
    load_areas::Vector{<:Real}      # Area of strips [-³]
    load_volumes::Vector{<:Real}    # Volume of strips [-³]
    max_spans::Vector{<:Real}     # Spans of each cell [-]
    slab_depths::Vector{<:Real}   # Depth of each cell [-]
    plot_context::PlotContext      # Plot context
    load_dictionary::Dict{Any, Vector{Asap.AbstractLoad}}  # Dictionary comparing elements to loads
    trib_dictionary::Dict{Any, Any}  # Dictionary mapping elements to vector of coordinate tuples
    record_tributaries::Bool            # Whether to record lines
    slab_units::Symbol              # Unit of measurement for slab (:m, :mm, :in, :ft)
    raster_df::DataFrame            # DataFrame of raster points
    element_id_lookup_df::Dict{Tuple{Int,Int}, Element} # DataFrame of element id lookup
    i_holes::Vector{Vector{Int64}}  # Indices of holes in the model
    i_perimeter::Vector{Int64}      # Indices of perimeter in the model
    
    function SlabAnalysisParams(model::Asap.Model; 
                              slab_name::String="", 
                              slab_type::Symbol=:isotropic, 
                              load_type::Symbol=:determinate,
                              vector_1d::Vector{<:Real}=[1.0, 0.0],
                              slab_thickness::Real=0.0,
                              perp::Bool=false,
                              perp_vector_1d::Vector{<:Real}=[0.0, -1.0],
                              slab_sizer::Symbol=:cellular,
                              fix_param::Bool=true,
                              spacing::Float64=0.1,
                              area::Float64=-1.0,
                              areas::Vector{<:Real}=Float64[],
                              load_areas::Vector{<:Real}=Float64[],
                              load_volumes::Vector{<:Real}=Float64[],
                              max_spans::Vector{<:Real}=Float64[],
                              slab_depths::Vector{<:Real}=Float64[],
                              plot_analysis::Bool=false,
                              load_dictionary::Dict{Any, Vector{Asap.AbstractLoad}}=Dict{Any, Vector{Asap.AbstractLoad}}(),
                              trib_dictionary::Dict{Any, Any}=Dict{Any, Any}(),
                              record_tributaries::Bool=false,
                              slab_units::Symbol=:m,
                              raster_df::DataFrame=DataFrame(),
                              element_id_lookup_df::Dict{Tuple{Int,Int}, Element}=Dict{Tuple{Int,Int}, Element}(),
                              i_holes::Vector{Vector{Int64}}=Vector{Int64}[],
                              i_perimeter::Vector{Int64}=Int64[])

        @assert (slab_type in [:isotropic, :uniaxial, :orth_biaxial]) "Invalid slab type."
        @assert (slab_sizer in [:cellular, :uniform]) "Invalid slab sizing method."
        @assert (load_type in [:determinate, :indeterminate]) "Invalid load type."
        
        plot_context = PlotContext(plot_analysis, nothing, nothing)
        vector_1d = Float64.(vector_1d)
        element_id_lookup_df = get_element_id(model)

        new(model, slab_name, slab_type, load_type, vector_1d, slab_thickness, perp, perp_vector_1d, slab_sizer, fix_param, spacing, area, areas, 
            load_areas, load_volumes, max_spans, slab_depths, plot_context, load_dictionary, trib_dictionary, record_tributaries, slab_units, raster_df, element_id_lookup_df, i_holes, i_perimeter)
    end
end


"""
    SlabSizingParams

Parameters for beam sizing, including depth, unit, and various flags for sizing constraints.
"""
mutable struct SlabSizingParams
    model::Union{Asap.Model, Nothing}               # The model

    # input values
    live_load::Real              # Live load [load/area]
    superimposed_dead_load::Real  # Superimposed dead load [load/area]
    slab_dead_load::Real         # Slab dead load [load/area]
    façade_load::Real            # Façade load [load/area]
    live_factor::Real             # Live load factor [-]
    dead_factor::Real            # Dead load factor [-]
    beam_sizer::Symbol              # :discrete or :continuous
    max_depth::Real                  # Maximum allowable depth for the assembly
    beam_units::Symbol              # Unit of measurement for beams (:m, :mm, :in, :ft)

    # default input values
    max_assembly_depth::Bool         # Whether there is a maximum assembly depth
    deflection_limit::Bool            # Whether to apply deflection limits
    minimum_continuous::Bool         # Whether there is a minimum optimization value
    collinear::Union{Bool,Nothing}   # True if collinear members are sized together
    drawn::Bool                      # True if the model has been drawn
    element_ids::Vector{Int64}      # Element ids
    serviceability_lim::Real       # Divide length by limit to get serviceability
    catalog_discrete::Any           # Catalog for discrete sizing
    n_max_sections::Int            # Maximum number of sections per beam

    # calculated values
    area::Real                      # Area of the slab in beam units
    w::Real                      # Distributed load [load/area]
    self_weight::Vector{<:Real}    # Weight of each cell [load/area] -- depends on your input units
    max_beam_depth::Real             # Maximum allowable height for the beams
    M_maxs::Vector{<:Real}          # Maximum moments
    V_maxs::Vector{<:Real}          # Maximum shear forces
    x_maxs::Vector{<:Real}          # Maximum x-coordinates
    load_dictionary::Dict{Any, Vector{Asap.AbstractLoad}}  # Dictionary comparing elements to loads
    load_df::DataFrame             # DataFrame of load information
    
    # results
    minimizers::Vector{Vector}      # Minimum values
    minimums::Vector{<:Real}        # Minimum values
    ids::Vector{String}         # IDs   
    collinear_minimizers::Vector{Vector}  # Minimizers if collinear
    collinear_ids::Vector{String}  # Ids if collinear
    collinear_minimums::Vector{<:Real}  # Minimums if collinear

    # setup values
    verbose::Bool                    # Whether to print detailed output
    
    function SlabSizingParams(;
        model::Union{Asap.Model, Nothing}=nothing,
        
        # input values
        live_load::Real=0.0,              # Live load [load/area]
        superimposed_dead_load::Real=0.0,  # Superimposed dead load [load/area] 
        slab_dead_load::Real=0.0,         # Slab dead load [load/area]
        façade_load::Real=0.0,            # Façade load [load/area]
        live_factor::Real=1.0,             # Live load factor [-]
        dead_factor::Real=1.0,            # Dead load factor [-]
        beam_sizer::Symbol=:discrete,         # :discrete or :continuous
        max_depth::Real=0.0,                  # Maximum allowable depth for the assembly

        # default input values
        max_assembly_depth::Bool=true,        # Whether there is a maximum assembly depth
        deflection_limit::Bool=true,          # Whether to apply deflection limits
        minimum_continuous::Bool=true,        # Whether there is a minimum optimization value
        collinear::Union{Bool,Nothing}=false, # True if collinear members are sized together
        drawn::Bool=false,
        element_ids::Vector{Int64}=Int64[],
        beam_units::Symbol=:in,              # Unit of measurement for beams
        serviceability_lim::Real=360,    # Divide length by limit to get serviceability
        catalog_discrete::Any=allW_imperial(), # Catalog for discrete sizing
        n_max_sections::Int=0,              # Maximum number of sections per beam

        # calculated values
        area::Real=0.0,                      # Area of the slab in beam units
        w::Real=0.0,                      # Distributed load [load/area]
        self_weight::Vector{<:Real}=Float64[], # Weight of each cell [load/area]
        max_beam_depth::Real=0.0,             # Maximum allowable height for the beams
        M_maxs::Vector{<:Real}=Float64[],    # Maximum moments
        V_maxs::Vector{<:Real}=Float64[],    # Maximum shear forces
        x_maxs::Vector{<:Real}=Float64[],    # Maximum x-coordinates
        load_dictionary::Dict{Any, Vector{Asap.AbstractLoad}}=Dict{Any, Vector{Asap.AbstractLoad}}(), # Dictionary comparing elements to loads
        load_df::DataFrame=DataFrame(), # DataFrame of load information

        # results
        minimizers::Vector{Vector}=Vector[],  # Minimum values
        minimums::Vector{<:Real}=Float64[],  # Minimum values
        ids::Vector{String}=String[],         # IDs
        collinear_minimizers::Vector{Vector}=Vector[], # Minimizers if collinear
        collinear_ids::Vector{String}=String[],       # Ids if collinear
        collinear_minimums::Vector{<:Real}=Float64[], # Minimums if collinear

        # setup values
        verbose::Bool=false                   # Whether to print detailed output
    )
        @assert (beam_sizer in [:discrete, :continuous]) "Invalid beam sizing method."

        new(model, live_load, superimposed_dead_load, slab_dead_load, façade_load, live_factor, dead_factor, beam_sizer, max_depth,
            beam_units, max_assembly_depth, deflection_limit, minimum_continuous, collinear, drawn, element_ids,
            serviceability_lim, catalog_discrete, n_max_sections, area, w, self_weight, max_beam_depth, M_maxs, V_maxs,
            x_maxs, load_dictionary, load_df, minimizers, minimums, ids, collinear_minimizers, collinear_ids,
            collinear_minimums, verbose)
    end
end


"""
    SlabOptimResults

Results of slab optimization, including mass, carbon footprint, and force data.
"""
mutable struct SlabOptimResults <: AbstractOptParams
    slab_name::String               # Slab name
    slab_type::Symbol               # Slab type
    vector_1d::Vector{<:Real}      # Direction of uniaxial slab
    slab_sizer::Symbol              # Slab sizer
    beam_sizer::Symbol              # Beam sizer
    area::Float64                   # Area
    minimizers::Vector{Vector}      # Imperial
    minimums::Vector{<:Real}       # Imperial
    ids::Vector{String}             # Identifications
    collinear::Union{Bool,Nothing}  # True/false collinear
    max_depth::Real                 # Maximum depth [in]
    mass_beams::Float64             # Mass of beams [kg]
    norm_mass_beams::Float64        # Normalized mass of beams [kg/m^2]
    embodied_carbon_beams::Float64  # Embodied carbon of beams [kg]
    mass_slab::Float64              # Mass of slab [kg]
    norm_mass_slab::Float64         # Normalized mass of slab [kg/m^2]
    embodied_carbon_slab::Float64   # Embodied carbon of slab [kg]
    mass_rebar::Float64             # Mass of rebar [kg]
    norm_mass_rebar::Float64        # Normalized mass of rebar [kg/m^2]
    embodied_carbon_rebar::Float64   # Embodied carbon of rebar [kg]
    areas::Vector{<:Real}           # Areas [in^2]
    x::Vector{Vector}                 # X-coordinates [m]
    P::Vector{Vector}                 # Loads [kN]
    Px::Vector{Vector}                # X-coordinates associated with loads [m]
    My::Vector{Vector}              # Moments [kNm]
    Mn::Vector{<:Real}             # Nominal moments [kNm]
    Vy::Vector{Vector}              # Shear forces [kN]
    Vn::Vector{<:Real}             # Nominal shear forces [kN]
    Δ_local::Vector{Vector}         # Local displacements [m]
    Δ_global::Vector{Vector}        # Global displacements [m]
    sections::Vector{Any}           # Sections

    function SlabOptimResults(slab_name::String, slab_type::Symbol, vector_1d::Vector{<:Real}, slab_sizer::Symbol, beam_sizer::Symbol, area::Float64, minimizers::Vector{Vector}, minimums::Vector{<:Real}, ids::Vector{String}, collinear::Union{Bool,Nothing}, max_depth::Real, mass_beams::Float64, norm_mass_beams::Float64, embodied_carbon_beams::Float64, mass_slab::Float64, norm_mass_slab::Float64, embodied_carbon_slab::Float64, mass_rebar::Float64, norm_mass_rebar::Float64, embodied_carbon_rebar::Float64, areas::Vector{<:Real}, x::Vector{Vector}, P::Vector{Vector}, Px::Vector{Vector}, My::Vector{Vector}, Mn::Vector{<:Real}, Vy::Vector{Vector}, Vn::Vector{<:Real}, Δ_local::Vector{Vector}, Δ_global::Vector{Vector}, sections::Vector{Any})
        new(slab_name, slab_type, vector_1d, slab_sizer, beam_sizer, area, minimizers, minimums, ids, collinear, max_depth, mass_beams, norm_mass_beams, embodied_carbon_beams, mass_slab, norm_mass_slab, embodied_carbon_slab, mass_rebar, norm_mass_rebar, embodied_carbon_rebar, areas, x, P, Px, My, Mn, Vy, Vn, Δ_local, Δ_global, sections)
    end

    function SlabOptimResults(max_depth::Float64=0.0, collinear::Union{Bool,Nothing}=nothing)
        slab_name = ""
        slab_type = :uniaxial
        vector_1d = [1.0, 0.0]
        slab_sizer = :cellular
        beam_sizer = :continuous
        area = mass_beams = norm_mass_beams = embodied_carbon_beams = mass_slab = norm_mass_slab = embodied_carbon_slab = mass_rebar = norm_mass_rebar = embodied_carbon_rebar = 0.0
        minimizers = x = P = Px = My = Vy = Δ_local = Δ_global =sections = [[]]
        minimums = areas = Mn = Vn = [0.0]
        ids = [""]
        collinear = collinear
        max_depth = max_depth

        new(slab_name, slab_type, vector_1d, slab_sizer, beam_sizer, area, minimizers, minimums, ids, collinear, max_depth, mass_beams, norm_mass_beams, embodied_carbon_beams, mass_slab, norm_mass_slab, embodied_carbon_slab, mass_rebar, norm_mass_rebar, embodied_carbon_rebar, areas, x, P, Px, My, Mn, Vy, Vn, Δ_local, Δ_global, sections)
    end
end
