abstract type AbstractVariable end
abstract type IndependentVariable <: AbstractVariable end
abstract type AbstractOptParams end

"""
    PlotContext

Context for plotting.

# Fields
"""
mutable struct PlotContext
    plot::Bool   # Plotting flag
    fig::Any     # CairoMakie Figure (or nothing when headless / no plot)
    ax::Any      # CairoMakie Axis (or nothing when headless / no plot)

    function PlotContext(plot::Bool, fig=nothing, ax=nothing)
        new(plot, fig, ax)
    end
end

"""
    SlabAnalysisParams

Parameters for slab analysis, including pre-sizing and post-sizing attributes.

Use `slab_depth_minimum` (same units as `slab_units`) to clamp slab depths **after**
the span-based thickness (or fixed `slab_thickness`) is computed, before loads are
assembled. The default ``0`` is stored as ``0.001`` m (in your `slab_units`) so depths
never sit below that numerical floor; larger requests are kept and also cannot fall
below ``0.001`` m.
"""
mutable struct SlabAnalysisParams <: AbstractOptParams
    model::Asap.Model               # The model
    slab_name::String               # Slab name
    slab_type::Symbol              # :uniaxial or :isotropic
    load_type::Symbol              # :determinate or :indeterminate
    vector_1d::Vector{<:Real}      # Direction of uniaxial slab
    slab_thickness::Real          # Slab thickness [slab length units, typically m]
    """Minimum slab depth after span-based thickness. Same units as `slab_thickness`.
    `0` means use the default floor (`0.001` m in `slab_units`). Any other value is
    raised to at least that same floor."""
    slab_depth_minimum::Real
    perp::Bool                    # Whether to use the perpendicular direction
    perp_vector_1d::Vector{<:Real} # Perpendicular direction of orth_biaxial slab
    slab_sizer::Symbol             # :cellular or :uniform
    fix_param::Bool                # Adjust for stiffness of elements or stay at 0.5?
    spacing::Float64               # Density of sampling [-]
    area::Float64                  # Sum area of slab [-²]
    areas::Vector{<:Real}         # Areas of each cell [-²]
    load_areas::Vector{<:Real}      # Area of strips [-²]
    load_volumes::Vector{<:Real}    # Volume of strips [-³]
    load_widths::Vector{<:Real}     # Perpendicular tributary width at each strip sample point [-]
    max_spans::Vector{<:Real}     # Spans of each cell [-]
    slab_depths::Vector{<:Real}   # Depth of each cell [-]
    plot_context::PlotContext      # Plot context
    load_dictionary::Dict{Tuple{Int,Int}, Vector{Asap.AbstractLoad}}  # Dictionary comparing elements to loads
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
                              slab_depth_minimum::Real=0.0,
                              perp::Bool=false,
                              perp_vector_1d::Vector{<:Real}=[0.0, -1.0],
                              slab_sizer::Symbol=:cellular,
                              fix_param::Bool=true,
                              spacing::Float64=0.1,
                              area::Float64=-1.0,
                              areas::Vector{<:Real}=Float64[],
                              load_areas::Vector{<:Real}=Float64[],
                              load_volumes::Vector{<:Real}=Float64[],
                              load_widths::Vector{<:Real}=Float64[],
                              max_spans::Vector{<:Real}=Float64[],
                              slab_depths::Vector{<:Real}=Float64[],
                              plot_analysis::Bool=false,
                              load_dictionary::Dict{Tuple{Int,Int}, Vector{Asap.AbstractLoad}}=Dict{Tuple{Int,Int}, Vector{Asap.AbstractLoad}}(),
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
        @assert slab_depth_minimum >= 0 "slab_depth_minimum must be ≥ 0."
        _floor_in_slab_units = 0.001 / convert_to_m[slab_units]
        slab_depth_minimum = max(slab_depth_minimum, _floor_in_slab_units)

        plot_context = PlotContext(plot_analysis, nothing, nothing)
        vector_1d = Float64.(vector_1d)
        element_id_lookup_df = get_element_id(model)

        new(model, slab_name, slab_type, load_type, vector_1d, slab_thickness, slab_depth_minimum, perp, perp_vector_1d, slab_sizer, fix_param, spacing, area, areas, 
            load_areas, load_volumes, load_widths, max_spans, slab_depths, plot_context, load_dictionary, trib_dictionary, record_tributaries, slab_units, raster_df, element_id_lookup_df, i_holes, i_perimeter)
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
    façade_load::Real            # Façade line load [load/length], typically kip/in
    live_factor::Real             # Live load factor [-]
    dead_factor::Real            # Dead load factor [-]
    beam_sizer::Symbol              # :discrete or :continuous
    nlp_solver::Symbol              # :Ipopt, :MMA (NLopt LD_MMA), :SLSQP, :CCSAQ, :COBYLA
    max_depth::Real                  # Maximum allowable depth for the assembly
    beam_units::Symbol              # Unit of measurement for beams (:m, :mm, :in, :ft)

    # default input values
    max_assembly_depth::Bool         # Whether there is a maximum assembly depth
    deflection_limit::Bool            # Whether to apply deflection limits
    staged_deflection_limit::Bool     # Whether to apply staged L/360 + L/240 checks when staged loads are available
    minimum_continuous::Bool         # Whether there is a minimum optimization value
    collinear::Union{Bool,Nothing}   # True if collinear members are sized together
    drawn::Bool                      # True if the model has been drawn
    element_ids::Vector{Int64}      # Element ids
    serviceability_lim::Real       # Divide length by limit to get serviceability
    catalog_discrete               # Catalog for discrete sizing
    n_max_sections::Int            # Maximum number of sections per beam

    # calculated values
    area::Real                      # Area of the slab in beam units
    w::Real                      # Distributed load [load/area]
    self_weight::Vector{<:Real}    # Weight of each cell [load/area] -- depends on your input units
    max_beam_depth::Real             # Maximum allowable height for the beams
    M_maxs::Vector{<:Real}          # Maximum moments
    V_maxs::Vector{<:Real}          # Maximum shear forces
    x_maxs::Vector{<:Real}          # Maximum x-coordinates
    load_dictionary::Dict{Tuple{Int,Int}, Vector{Asap.AbstractLoad}}  # Dictionary comparing elements to loads
    load_df::DataFrame             # DataFrame of load information
    
    # results
    minimizers::Vector{Vector{Float64}}      # Minimum values
    minimums::Vector{Float64}        # Minimum values
    ids::Vector{String}         # IDs   
    collinear_minimizers::Vector{Vector{Float64}}  # Minimizers if collinear
    collinear_ids::Vector{String}  # Ids if collinear
    collinear_minimums::Vector{Float64}  # Minimums if collinear
    collinear_groups::Vector{Int}  # Group ID per beam element (populated during sizing)

    # composite action
    composite_action::Bool           # Use composite stiffness for deflection checks
    E_c::Real                        # Concrete elastic modulus (ksi) for modular ratio
    slab_depth_in::Real              # Slab depth in beam units (in), set during sizing
    deflection_reduction_factor::Real # Divide computed deflection by this before checking limit (e.g. 2.5 for bare steel w/ slab)
    i_perimeter::Set{Int}            # Beam element indices on the slab perimeter (edge beams)

    # staged deflection Ix floors (populated by outer convergence loop)
    min_Ix_comp::Dict{Int,Float64}   # Per-beam composite Ix lower bound from staged verification
    min_Ix_bare::Dict{Int,Float64}   # Per-beam bare-steel Ix lower bound from staged verification

    # global deflection sanity check
    max_bay_span::Real               # Largest column bay span in beam units (in), for global δ warning

    # MIP baseline (populated when NLP runs MIP internally as a warm-start)
    mip_result::Union{SlabSizingParams, Nothing}

    # staged deflection convergence (populated by outer loop in optimal_beamsizer)
    staged_converged::Bool           # true if staged deflection loop converged (or was not applicable)
    staged_n_violations::Int         # number of beams still violating limits at exit (0 = converged)

    # setup values
    verbose::Bool                    # Whether to print detailed output

    # material scenario (used by load generation + postprocessing for density/ECC)
    concrete_material::ConcreteMaterial

    # factored beam self-weight as `LineLoad`s on FE model (strength / columns); empty until sizing
    beam_sw_line_loads::Vector{Asap.AbstractLoad}
    
    function SlabSizingParams(;
        model::Union{Asap.Model, Nothing}=nothing,
        
        # input values
        live_load::Real=0.0,              # Live load [load/area]
        superimposed_dead_load::Real=0.0,  # Superimposed dead load [load/area] 
        slab_dead_load::Real=0.0,         # Slab dead load [load/area]
        façade_load::Real=0.0,            # Façade line load [load/length], typically kip/in
        live_factor::Real=1.0,             # Live load factor [-]
        dead_factor::Real=1.0,            # Dead load factor [-]
        beam_sizer::Symbol=:discrete,         # :discrete or :continuous
        nlp_solver::Symbol=:Ipopt,             # :Ipopt, :MMA (NLopt LD_MMA), :SLSQP, :CCSAQ, :COBYLA
        max_depth::Real=0.0,                  # Maximum allowable depth for the assembly

        # default input values
        max_assembly_depth::Bool=true,        # Whether there is a maximum assembly depth
        deflection_limit::Bool=true,          # Whether to apply deflection limits
        staged_deflection_limit::Bool=true,   # Whether to apply staged L/360 + L/240 checks when staged loads are available
        minimum_continuous::Bool=true,        # Whether there is a minimum optimization value
        collinear::Union{Bool,Nothing}=false, # True if collinear members are sized together
        drawn::Bool=false,
        element_ids::Vector{Int64}=Int64[],
        beam_units::Symbol=:in,              # Unit of measurement for beams
        serviceability_lim::Real=360,    # Divide length by limit to get serviceability
        catalog_discrete=allW_imperial(), # Catalog for discrete sizing
        n_max_sections::Int=0,              # Maximum number of sections per beam

        # calculated values
        area::Real=0.0,                      # Area of the slab in beam units
        w::Real=0.0,                      # Distributed load [load/area]
        self_weight::Vector{<:Real}=Float64[], # Weight of each cell [load/area]
        max_beam_depth::Real=0.0,             # Maximum allowable height for the beams
        M_maxs::Vector{<:Real}=Float64[],    # Maximum moments
        V_maxs::Vector{<:Real}=Float64[],    # Maximum shear forces
        x_maxs::Vector{<:Real}=Float64[],    # Maximum x-coordinates
        load_dictionary::Dict{Tuple{Int,Int}, Vector{Asap.AbstractLoad}}=Dict{Tuple{Int,Int}, Vector{Asap.AbstractLoad}}(), # Dictionary comparing elements to loads
        load_df::DataFrame=DataFrame(), # DataFrame of load information

        # results
        minimizers::Vector{Vector{Float64}}=Vector{Float64}[],  # Minimum values
        minimums::Vector{Float64}=Float64[],  # Minimum values
        ids::Vector{String}=String[],         # IDs
        collinear_minimizers::Vector{Vector{Float64}}=Vector{Float64}[], # Minimizers if collinear
        collinear_ids::Vector{String}=String[],       # Ids if collinear
        collinear_minimums::Vector{Float64}=Float64[], # Minimums if collinear
        collinear_groups::Vector{Int}=Int[],           # Group ID per beam element

        # composite action
        composite_action::Bool=false,         # Use composite stiffness for deflection
        E_c::Real=57.0 * sqrt(4000.0),          # Concrete E_c (ksi), = 57√(f'c_psi) for f'c = 4 ksi
        slab_depth_in::Real=0.0,              # Slab depth (in), populated during sizing
        deflection_reduction_factor::Real=1.0, # Divide computed deflection by this (e.g. 2.5 for bare steel w/ slab)
        i_perimeter::Set{Int}=Set{Int}(),     # Perimeter beam indices, populated during sizing

        # staged deflection Ix floors (populated by outer convergence loop)
        min_Ix_comp::Dict{Int,Float64}=Dict{Int,Float64}(),
        min_Ix_bare::Dict{Int,Float64}=Dict{Int,Float64}(),

        # global deflection sanity check
        max_bay_span::Real=0.0,              # Largest column bay span in beam units (in)

        # MIP baseline (auto-populated when NLP runs MIP as warm-start)
        mip_result::Union{SlabSizingParams, Nothing}=nothing,

        # staged deflection convergence
        staged_converged::Bool=true,
        staged_n_violations::Int=0,

        # setup values
        verbose::Bool=false,                  # Whether to print detailed output

        # material scenario
        concrete_material::ConcreteMaterial=DEFAULT_CONCRETE,

        beam_sw_line_loads::Vector{Asap.AbstractLoad}=Asap.AbstractLoad[],
    )
        @assert (beam_sizer in [:discrete, :continuous]) "Invalid beam sizing method."
        @assert deflection_reduction_factor > 0 "deflection_reduction_factor must be > 0."

        # High-fidelity policy: when composite stiffness is explicitly modeled, do not
        # also reduce deflection by an empirical factor (avoids double-counting stiffness).
        drf = composite_action ? 1.0 : deflection_reduction_factor

        new(model, live_load, superimposed_dead_load, slab_dead_load, façade_load, live_factor, dead_factor, beam_sizer, nlp_solver, max_depth,
            beam_units, max_assembly_depth, deflection_limit, staged_deflection_limit, minimum_continuous, collinear, drawn, element_ids,
            serviceability_lim, catalog_discrete, n_max_sections, area, w, self_weight, max_beam_depth, M_maxs, V_maxs,
            x_maxs, load_dictionary, load_df, minimizers, minimums, ids, collinear_minimizers, collinear_ids,
            collinear_minimums, collinear_groups, composite_action, E_c, slab_depth_in, drf, i_perimeter,
            min_Ix_comp, min_Ix_bare, max_bay_span, mip_result,
            staged_converged, staged_n_violations, verbose, concrete_material, beam_sw_line_loads)
    end
end


"""
    SlabOptimResults

Results of slab optimization, including mass, embodied carbon, internal forces,
column sizing, and staged deflection analysis.

Uses `@kwdef` so fields can be set by name — adding new fields only requires
a default value, not updating every call site.

# Staged deflection fields

For unshored composite construction, deflections are decomposed by load stage:

| Field           | Description                                          |
|-----------------|------------------------------------------------------|
| `δ_slab_dead`   | Slab DL deflection on bare steel Ix [in]           |
| `δ_beam_dead`   | Beam self-weight deflection on bare steel [in]     |
| `δ_sdl`         | SDL deflection on composite Ix [in]                |
| `δ_live`        | Live load deflection on composite Ix [in]          |
| `δ_total`       | Superposition of all stages [in]                   |
| `Δ_limit_live`  | L/360 limit per beam [in]                          |
| `Δ_limit_total` | L/240 limit per beam [in]                          |
| `δ_live_ok`     | true if δ_live ≤ L/360                              |
| `δ_total_ok`    | true if δ_total ≤ L/240                             |

Sizer exit state (from `optimal_beamsizer` outer staged-deflection loop):

| Field                 | Description                                        |
|-----------------------|----------------------------------------------------|
| `composite_action`    | Whether composite stiffness was used in sizing     |
| `staged_converged`    | `true` if staged Ix loop converged (or N/A)        |
| `staged_n_violations` | Beams still violating staged limits at exit        |
| `staged_ok`          | Convenience summary: staged loop exited cleanly    |
"""
@kwdef mutable struct SlabOptimResults <: AbstractOptParams
    slab_name::String                          = ""
    slab_type::Symbol                          = :uniaxial
    vector_1d::Vector{Float64}                 = [1.0, 0.0]
    slab_sizer::Symbol                         = :cellular
    beam_sizer::Symbol                         = :continuous
    area::Float64                              = 0.0
    minimizers::Vector{Vector{Float64}}        = Vector{Float64}[Float64[]]
    minimums::Vector{Float64}                  = [0.0]
    ids::Vector{String}                        = [""]
    collinear::Union{Bool,Nothing}             = nothing
    max_depth::Float64                         = 0.0
    mass_beams::Float64                        = 0.0
    norm_mass_beams::Float64                   = 0.0
    embodied_carbon_beams::Float64             = 0.0
    mass_slab::Float64                         = 0.0
    norm_mass_slab::Float64                    = 0.0
    embodied_carbon_slab::Float64              = 0.0
    mass_rebar::Float64                        = 0.0
    norm_mass_rebar::Float64                   = 0.0
    embodied_carbon_rebar::Float64             = 0.0
    areas::Vector{Float64}                     = [0.0]
    x::Vector{Vector{Float64}}                 = Vector{Float64}[Float64[]]
    P::Vector{Vector{Float64}}                 = Vector{Float64}[Float64[]]
    Px::Vector{Vector{Float64}}                = Vector{Float64}[Float64[]]
    My::Vector{Vector{Float64}}                = Vector{Float64}[Float64[]]
    Mn::Vector{Float64}                        = [0.0]
    Vy::Vector{Vector{Float64}}                = Vector{Float64}[Float64[]]
    Vn::Vector{Float64}                        = [0.0]
    Δ_local::Vector{Vector{Float64}}           = Vector{Float64}[Float64[]]
    Δ_global::Vector{Vector{Float64}}          = Vector{Float64}[Float64[]]
    sections::Vector{String}                   = [""]

    # --- Column results ---
    col_sections::Vector{String}               = String[]
    col_Pu::Vector{Float64}                    = Float64[]
    col_ϕPn::Vector{Float64}                   = Float64[]
    col_util::Vector{Float64}                  = Float64[]
    mass_columns::Float64                      = 0.0
    norm_mass_columns::Float64                 = 0.0
    embodied_carbon_columns::Float64           = 0.0

    # --- Fireproofing results ---
    mass_fireproofing::Float64                 = 0.0
    norm_mass_fireproofing::Float64            = 0.0
    embodied_carbon_fireproofing::Float64      = 0.0

    # --- Staged deflection results (per beam) ---
    δ_slab_dead::Vector{Float64}               = Float64[]
    δ_beam_dead::Vector{Float64}               = Float64[]
    δ_sdl::Vector{Float64}                     = Float64[]
    δ_live::Vector{Float64}                    = Float64[]
    δ_total::Vector{Float64}                   = Float64[]
    Δ_limit_live::Vector{Float64}              = Float64[]
    Δ_limit_total::Vector{Float64}             = Float64[]
    δ_live_ok::Vector{Bool}                    = Bool[]
    δ_total_ok::Vector{Bool}                   = Bool[]

    # --- Global deflection sanity check ---
    max_δ_total::Float64                       = 0.0
    max_bay_span::Float64                      = 0.0
    global_δ_ok::Bool                          = true

    # --- Aggregate utilization ---
    max_util_M::Float64                        = 0.0
    max_util_V::Float64                        = 0.0
    max_col_util::Float64                      = 0.0
    n_L360_fail::Int                           = 0
    i_L360_fail::Vector{Int}                   = Int[]
    n_L240_fail::Int                           = 0
    i_L240_fail::Vector{Int}                   = Int[]

    # --- Optimizer / serviceability flags (from `SlabSizingParams`) ---
    nlp_solver::String                         = ""   # `"MIP"` for discrete; else e.g. `"MMA"`, `"Ipopt"`
    deflection_limit::Bool                     = true
    staged_deflection_limit::Bool              = true

    # --- Sizer staged-deflection loop (see `optimal_beamsizer`) ---
    composite_action::Bool                     = false
    staged_converged::Bool                     = true
    staged_n_violations::Int                   = 0
    staged_ok::Bool                            = true

    # --- Config versioning (for resume-logic staleness detection) ---
    config_hash::String                        = ""

    # --- Run diagnostics (for CSV reporting / filtering) ---
    geometry_file::String                      = ""
    span_ok::Bool                              = false
    result_ok::Bool                            = false
    strength_ok::Bool                          = false
    serviceability_ok::Bool                    = false
    column_ok::Bool                            = false
    solver_status::String                      = ""
    diagnostic_flags::String                   = ""
    diagnostic_messages::String                = ""
end
