
"""
    reset_SlabAnalysisParams(self::SlabAnalysisParams, geometry_dict::Dict; plot::Bool=false, drawn::Bool=false, sections::Vector=[])

Resets a slab analysis object with new geometry while preserving original parameters.

# Arguments
- `self::SlabAnalysisParams`: The original slab analysis parameters
- `geometry_dict::Dict`: Dictionary containing new geometry data
- `plot::Bool=false`: Whether to plot the geometry
- `drawn::Bool=false`: Whether the geometry has been drawn
- `sections::Vector=[]`: Initial sections for the beams

# Returns
- A new `SlabAnalysisParams` object with reset geometry but preserved parameters

"""
function reset_SlabAnalysisParams(self::SlabAnalysisParams, model::Asap.Model; plot::Bool=false, drawn::Bool=false, sections::Vector=[])
    return SlabAnalysisParams(
        model,
        slab_name=self.slab_name,
        slab_type=self.slab_type,
        vector_1d=self.vector_1d,
        slab_sizer=self.slab_sizer,
        fix_param=self.fix_param,
        spacing=self.spacing,
        area=self.area,
        areas=self.areas,
        load_areas=self.load_areas,
        load_volumes=self.load_volumes,
        max_spans=self.max_spans,
        slab_depths=self.slab_depths,
        plot_analysis=self.plot_context.plot,
        load_dictionary=self.load_dictionary,
        trib_dictionary=self.trib_dictionary,
        slab_units=self.slab_units,
        i_holes=self.i_holes,
        i_perimeter=self.i_perimeter,
    )
end

"""
    reset_SlabSizingParams(self::SlabSizingParams)

Resets a slab sizing object with new parameters while preserving original geometry.
"""
function reset_SlabSizingParams(self::SlabSizingParams)
    return SlabSizingParams(
        model=self.model,
        live_load=self.live_load,
        superimposed_dead_load=self.superimposed_dead_load,
        live_factor=self.live_factor,
        dead_factor=self.dead_factor,
        beam_sizer=self.beam_sizer,
        max_depth=self.max_depth,
        max_assembly_depth=self.max_assembly_depth,
        deflection_limit=self.deflection_limit,
        minimum_continuous=self.minimum_continuous,
        collinear=self.collinear,
        beam_units=self.beam_units,
        serviceability_lim=self.serviceability_lim,
        catalog_discrete=self.catalog_discrete,
        w=self.w,
        self_weight=self.self_weight,
        max_beam_depth=self.max_beam_depth,
        M_maxs=Float64[],
        V_maxs=Float64[],
        x_maxs=Float64[],
        load_dictionary=self.load_dictionary,
        load_df=self.load_df,
        minimizers=Vector{Vector}(),
        minimums=Float64[],
        ids=String[],
        collinear_minimizers=Vector{Vector}(),
        collinear_ids=String[],
        collinear_minimums=Float64[],
        verbose=self.verbose
    )
end