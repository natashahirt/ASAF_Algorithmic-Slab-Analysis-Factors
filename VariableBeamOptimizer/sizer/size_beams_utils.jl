"""
    get_deflection_constraint(beam_params, beam_length, params; beam_element, element_loads)

Single serviceability inequality constraint for the beam NLP.

# Composite branch (`params.composite_action == true`)

Applies the **total unfactored service load** (slab DL + SDL + LL + beam
self-weight) to the **transformed composite section** and checks

    δ_composite = δ_bare · (I_bare / I_composite)  ≤  L / `serviceability_lim`

where

- `δ_bare` is the FE deflection of the bare steel beam under the total
  unfactored load (the FE adds self-weight via `dead_load = steel.ρ`), and
- `I_composite` is the effective transformed-section moment of inertia per
  AISC 360 §I3.1a (`get_I_composite_effective`), with concrete transformed
  by the modular ratio `n = E_s / E_c` and the effective width bounded by
  the `L/8` per-side cap.

This is a deliberate single-limit preliminary-sizing formulation: no
staged (L/360-live + L/240-total) split, no construction-stage bookkeeping,
no deflection-reduction factor. The rationale and citations are documented
in the accompanying paper (see `size_beams.jl::optimal_beamsizer`).

# Bare-steel branch (`params.composite_action == false`)

Bare-steel Ix with analytical self-weight contribution, divided by
`params.deflection_reduction_factor` (typically 1.0 for pure bare-steel
design; kept for Nov 2024 compatibility where a composite-stiffness
multiplier was used instead of a transformed section):

    (δ_bare + δ_SW) / DRF  ≤  L / `serviceability_lim`

# Returns

Closure `vars -> Vector{Float64}` of length 1 with a single residual
`δ - δ_max` (feasible when `≤ 0`). Differentiable in `vars` — the
composite `Ix` is recomputed each call.
"""
function get_deflection_constraint(beam_params::FrameOptParams, beam_length::Real, params::SlabSizingParams;
                                   beam_element::Union{Element,Nothing}=nothing,
                                   element_loads::Union{Vector,Nothing}=nothing)

    use_composite = params.composite_action && !isnothing(beam_element) && !isnothing(element_loads)

    composite_ctx = nothing
    if use_composite
        loadid_index = _build_loadid_index(params)
        load_positions   = Float64[]
        load_trib_widths = Float64[]
        for ld in element_loads
            if hasproperty(ld, :loadID)
                row = get(loadid_index, getproperty(ld, :loadID), nothing)
                if !isnothing(row)
                    push!(load_positions,   ld.position)
                    push!(load_trib_widths, params.load_df[row, :trib_width])
                end
            end
        end
        beam_idx = findfirst(el -> el === beam_element, params.model.elements[:beam])
        is_perim = !isnothing(beam_idx) && beam_idx in params.i_perimeter
        composite_ctx = (L_beam       = beam_element.length,
                         positions    = load_positions,
                         widths       = load_trib_widths,
                         t_slab       = params.slab_depth_in,
                         E_s          = steel_ksi.E,
                         E_c          = params.E_c,
                         is_perimeter = is_perim)
    end

    ρ_steel_kip = steel_ksi.ρ

    # Effective deflection limit after applying the reconciliation-loop
    # tightening factor. `serviceability_tighten ≥ 1` shrinks the allowable
    # δ so the sizer targets e.g. L/(360·1.15) internally while the paper
    # still reports L/360; set to 1.0 for a plain single-limit regime.
    δ_max = beam_length / (params.serviceability_lim * params.serviceability_tighten)

    # Composite stiffness knockdown (AISC partial shear connection / slip /
    # creep). Applied to `Ix_comp` in the sizer; `postprocess_slab` applies
    # the same factor so the two paths stay consistent.
    pcf = params.partial_composite_factor

    function single_section_deflection(vars::Vector)
        # FE deflection under total unfactored point loads + beam self-weight
        # (dead_load defaults to steel density when material is passed).
        δ_local = Zygote.ignore() do
            get_element_deflection(beam_params, vars, material=steel_ksi)
        end
        δ_bare = maximum(abs.(δ_local))

        h, w, tw, tf = vars[1], vars[2], vars[3], vars[4]
        Ix_bare = Ix_I_symm(h, w, tw, tf)

        if use_composite
            Ix_comp_full = get_I_composite_effective(h, w, tw, tf,
                         composite_ctx.t_slab, composite_ctx.E_s, composite_ctx.E_c,
                         composite_ctx.L_beam,
                         composite_ctx.positions, composite_ctx.widths;
                         is_perimeter=composite_ctx.is_perimeter)
            # Partial-composite knockdown: Ix_eff = Ix_bare + pcf · (Ix_comp − Ix_bare).
            # At pcf = 1 this is the full transformed Ix; at pcf = 0 it collapses
            # to bare steel. Preserves Ix_eff ≥ Ix_bare under any knockdown,
            # which the naive `pcf · Ix_comp` form does not.
            Ix_eff      = Ix_bare + pcf * (Ix_comp_full - Ix_bare)
            δ_composite = δ_bare * Ix_bare / Ix_eff
            return [δ_composite - δ_max]
        end

        # Bare-steel fallback: divide by DRF (kept for Nov 2024-style designs
        # that use a stiffness multiplier instead of a transformed section).
        return [δ_bare / params.deflection_reduction_factor - δ_max]
    end

    return single_section_deflection
end

"""
    check_uniqueness!(M_max, V_max, x_max, M_maxs, V_maxs, x_maxs, minimizers, minimums, ids)

Checks if the current section is unique and updates the lists if not.
"""
function check_uniqueness!(params::SlabSizingParams, M_max, V_max, x_max)
    unique = true
    if isempty(params.minimizers)
        return params, unique
    else
        for j in 1:lastindex(params.M_maxs) # check uniqueness within 1%
            if abs((params.M_maxs[j] - M_max)/M_max) < 0.01 && 
               abs((params.V_maxs[j] - V_max)/V_max) < 0.01 && 
               abs((params.x_maxs[j] - x_max)/x_max) < 0.01
                unique = false
                println("Not unique")
                push!(params.minimizers, copy(params.minimizers[j]))
                push!(params.minimums, copy(params.minimums[j]))
                push!(params.ids, params.ids[j])
                push!(params.M_maxs, M_max)
                push!(params.V_maxs, V_max)
                push!(params.x_maxs, x_max)
                return params, unique
            end
        end
        return params, unique
    end
end


"""
    initialize_slab_results(analysis_params::SlabAnalysisParams, sizing_params::SlabSizingParams)

Initializes the slab optimization results for both discrete and continuous beam sizing methods.

# Arguments
- `analysis_params::SlabAnalysisParams`: Parameters for slab analysis including slab name, type, etc.
- `sizing_params::SlabSizingParams`: Parameters for beam sizing including maximum depth constraints.

# Returns
- `slab_results_discrete_noncollinear`: Results for discrete sizing without collinear element grouping
- `slab_results_discrete_collinear`: Results for discrete sizing with collinear element grouping  
- `slab_results_continuous_noncollinear`: Results for continuous sizing without collinear element grouping
- `slab_results_continuous_collinear`: Results for continuous sizing with collinear element grouping
"""
function initialize_slab_results(analysis_params::SlabAnalysisParams, sizing_params::SlabSizingParams)

    # Initialize results
    slab_results_discrete_noncollinear = SlabOptimResults()
    slab_results_discrete_collinear = SlabOptimResults()
    slab_results_continuous_noncollinear = SlabOptimResults()
    slab_results_continuous_collinear = SlabOptimResults()

    # clean slab result
    for slab_result in [slab_results_discrete_noncollinear, slab_results_discrete_collinear, slab_results_continuous_noncollinear, slab_results_continuous_collinear]
        slab_result.slab_name = analysis_params.slab_name
        slab_result.slab_type = analysis_params.slab_type
        slab_result.vector_1d = analysis_params.vector_1d
        slab_result.slab_sizer = analysis_params.slab_sizer
        slab_result.max_depth = sizing_params.max_depth
    end

    slab_results_discrete_noncollinear.collinear = false
    slab_results_discrete_collinear.collinear = true
    slab_results_continuous_noncollinear.collinear = false
    slab_results_continuous_collinear.collinear = true

    slab_results_discrete_noncollinear.beam_sizer = :discrete
    slab_results_discrete_collinear.beam_sizer = :discrete
    slab_results_continuous_noncollinear.beam_sizer = :continuous
    slab_results_continuous_collinear.beam_sizer = :continuous

    return slab_results_discrete_noncollinear, slab_results_discrete_collinear, slab_results_continuous_noncollinear, slab_results_continuous_collinear

end

"""
    get_collinear_groups(beam_elements) -> Vector{Int}

Compute collinear group IDs for a vector of beam elements using union-find.
Adjacent beams that share a node and pass `check_collinearity` are placed in the
same group. Returns a vector of integer group IDs (1-indexed, contiguous).
"""
function get_collinear_groups(beam_elements)
    n = length(beam_elements)
    n == 0 && return Int[]

    # Build node → incident beam indices map and mark "break" nodes.
    # Break nodes terminate collinear grouping, e.g. supports/walls and beam junctions.
    node_to_beams = Dict{Any, Vector{Int}}()
    for i in 1:n
        ns = beam_elements[i].nodeStart.nodeID
        ne = beam_elements[i].nodeEnd.nodeID
        push!(get!(node_to_beams, ns, Int[]), i)
        push!(get!(node_to_beams, ne, Int[]), i)
    end

    break_nodes = Set{Any}()
    for (nid, incident) in node_to_beams
        # Topology branch points should split groups.
        if length(incident) != 2
            push!(break_nodes, nid)
            continue
        end
        # Explicit support-like IDs should split groups.
        if nid isa Symbol
            s = String(nid)
            if occursin("wall", s) || occursin("fixed", s) || occursin("column", s) || occursin("col", s)
                push!(break_nodes, nid)
            end
        end
    end

    parent = collect(1:n)
    function uf_find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    function uf_union(a, b)
        ra, rb = uf_find(a), uf_find(b)
        if ra != rb
            parent[rb] = ra
        end
    end

    for i in 1:n, j in (i+1):n
        ids_i = Asap.nodeids(beam_elements[i])
        ids_j = Asap.nodeids(beam_elements[j])
        shared = intersect(ids_i, ids_j)
        if !isempty(shared) && check_collinearity(beam_elements[i], beam_elements[j])
            # Do not group across support/junction/column break nodes.
            if any(nid -> nid in break_nodes, shared)
                continue
            end
            uf_union(i, j)
        end
    end

    roots = [uf_find(i) for i in 1:n]
    uid = sort(unique(roots))
    root_to_gid = Dict(r => g for (g, r) in enumerate(uid))
    return [root_to_gid[r] for r in roots]
end

"""
    collect_collinear_elements(self::SlabAnalysisParams, params::SlabSizingParams)

Ensures all beams in a collinear group share the same (heaviest) section.
Uses `get_collinear_groups` for grouping, then upgrades every member to the
section with the largest cross-sectional area.

# Returns
- `(minimizers, ids, minimums)` with uniform sections within each group.
"""
function collect_collinear_elements(self::SlabAnalysisParams, params::SlabSizingParams)
    beam_elements = self.model.elements[:beam]

    minimizers = deepcopy(params.minimizers)
    ids = deepcopy(params.ids)
    minimums = deepcopy(params.minimums)

    collinear_groups = get_collinear_groups(beam_elements)

    group_to_indices = Dict{Int, Vector{Int}}()
    for (idx, gid) in enumerate(collinear_groups)
        push!(get!(group_to_indices, gid, Int[]), idx)
    end

    for (_, element_idx) in group_to_indices
        length(element_idx) <= 1 && continue

        collinear_areas = [I_symm(vars...).A for vars in minimizers[element_idx]]
        _, max_idx = findmax(collinear_areas)

        for idx in element_idx
            minimizers[idx] = minimizers[element_idx[max_idx]]
            ids[idx] = ids[element_idx[max_idx]]
            minimums[idx] = minimums[element_idx[max_idx]]
        end
    end

    return minimizers, ids, minimums
end

"""
    get_scaled_model(self::SlabAnalysisParams, params::SlabSizingParams, conversion_factor::Float64)

Creates a scaled model with converted units and loads for analysis.

# Arguments
- `self::SlabAnalysisParams`: The slab analysis parameters containing the original model
- `params::SlabSizingParams`: Parameters for slab sizing including load factors and loads
- `conversion_factor::Float64`: Factor to convert from metric to imperial units (m to in)

# Returns
- Updates params.load_df with a DataFrame containing load information in imperial units
"""
function _normalize_section_to_inches(sec::Section)
    # Geometry imports often carry section properties in mm-based geometry units
    # even when E/G are already in ksi. Detect those large mm-scale properties and
    # convert them once so the scaled model is consistently kip-inch.
    looks_metric_geom = sec.A > 200 || sec.Ix > 1e6 || sec.Iy > 1e6 || sec.J > 1e6
    if !looks_metric_geom
        return sec
    end

    return Section(
        to_in(sec.A, power=2),  # mm² → in²
        sec.E,                  # ksi — already correct
        sec.G,                  # ksi — already correct
        to_in(sec.Ix, power=4), # mm⁴ → in⁴
        to_in(sec.Iy, power=4), # mm⁴ → in⁴
        to_in(sec.J, power=4)   # mm⁴ → in⁴
    )
end

function get_scaled_model(self::SlabAnalysisParams, params::SlabSizingParams, conversion_factor::Float64)
    # set up slab loads
    factored_w_applied = params.live_factor * params.live_load + params.dead_factor * params.superimposed_dead_load # ksi

    # Slab dead load basis (unfactored), then apply load factors consistently.
    ρ_conc_kipin3 = params.concrete_material.ρ_concrete_kipin3
    area_dead_load = false
    if params.slab_dead_load > 0.0
        unfactored_w_slab = params.slab_dead_load
        area_dead_load = true
    else
        if self.slab_type == :orth_biaxial
            unfactored_w_slab = 0.98 * ρ_conc_kipin3 + 0.02 * ρ_STEEL_KIPIN3
        else
            unfactored_w_slab = 0.99 * ρ_conc_kipin3 + 0.01 * ρ_STEEL_KIPIN3
        end
    end

    unfactored_w_applied = params.live_load + params.superimposed_dead_load
    factored_w_slab = params.dead_factor * unfactored_w_slab

    factored_w_façade = params.dead_factor * params.façade_load # kip/in
    unfactored_w_façade = params.façade_load # kip/in

    # Create a DataFrame to store load information.
    # `unfactored_w_live` and `unfactored_w_sdl` decompose `unfactored_w_applied`
    # for staged deflection calculations (bare steel DL vs composite SDL+LL).
    params.load_df = DataFrame(
        loadID = Int[],
        area = Float64[], # in²
        volume = Float64[], # in³
        depth = Float64[], # in
        trib_width = Float64[], # in — perpendicular tributary width at this strip
        factored_w_slab = Float64[], # kip
        unfactored_w_slab = Float64[], # kip
        factored_w_applied = Float64[], # kip
        unfactored_w_applied = Float64[], # kip
        unfactored_w_live = Float64[], # kip — live load component only
        unfactored_w_sdl = Float64[], # kip — superimposed dead load component only
    )

    # refactor the model without destroying the original
    factored_load_areas = self.load_areas * conversion_factor^2 # m² to in²
    factored_load_volumes = self.load_volumes * conversion_factor^3 # m³ to in³
    factored_load_widths = self.load_widths * conversion_factor # m to in
    factored_load_depths = factored_load_volumes ./ factored_load_areas # in³/in² -> in

    println("factored_load_depths: $(maximum(factored_load_depths)), $(minimum(factored_load_depths))")

    # start setting up scaled model
    new_nodes = Asap.Node[]
    new_loads = Asap.AbstractLoad[]
    new_elements = Asap.Element[]

    # apply appropriate slab loads
    for node in self.model.nodes
        new_position = node.position .* conversion_factor # m to in
        new_node = Asap.Node(new_position, node.dof)
        new_node.id = node.id
        push!(new_nodes, new_node)
    end

    for element in self.model.elements
        release = AsapToolkit.asap_release_symbol(element)
        sec = _normalize_section_to_inches(element.section)
        new_element = Asap.Element(new_nodes[element.nodeStart.nodeID], new_nodes[element.nodeEnd.nodeID], sec, release=release)
        new_element.id = element.id
        push!(new_elements, new_element)
    end

    element_to_index = Dict{Element, Int}()
    for (idx, el) in enumerate(self.model.elements)
        element_to_index[el] = idx
    end

    for (i, load) in enumerate(self.model.loads)

        if area_dead_load
            slab_load_value = factored_load_areas[i] * factored_w_slab
        else
            slab_load_value = factored_load_volumes[i] * factored_w_slab
        end

        applied_load_value = factored_load_areas[i] * factored_w_applied

        load_value = - (slab_load_value + applied_load_value)

        element_index = element_to_index[load.element]
        new_element = new_elements[element_index]
        new_load = PointLoad(new_element, load.position, [0,0,load_value])
        new_load.loadID = load.loadID
        push!(new_loads, new_load)

        # populate the dictionary — slab columns must store the same force
        # magnitude used for the point load (volume×density when density-based,
        # area×pressure when area-based)
        slab_factored_force   = area_dead_load ?
            factored_load_areas[i] * factored_w_slab :
            factored_load_volumes[i] * factored_w_slab
        slab_unfactored_force = area_dead_load ?
            factored_load_areas[i] * unfactored_w_slab :
            factored_load_volumes[i] * unfactored_w_slab

        push!(params.load_df, (new_load.loadID, 
                        factored_load_areas[i], 
                        factored_load_volumes[i], 
                        factored_load_depths[i], 
                        factored_load_widths[i],
                        slab_factored_force, 
                        slab_unfactored_force, 
                        factored_load_areas[i] * factored_w_applied, 
                        factored_load_areas[i] * unfactored_w_applied,
                        factored_load_areas[i] * params.live_load,
                        factored_load_areas[i] * params.superimposed_dead_load))
    end

    perimeter_set = Set(self.i_perimeter)
    for (i, element) in enumerate(self.model.elements[:beam])
        if i in perimeter_set
            # Attach perimeter facade loads to the scaled model elements, not the original model.
            scaled_idx = element_to_index[element]
            push!(new_loads, LineLoad(new_elements[scaled_idx], [0,0,-factored_w_façade]))
        end
    end

    # solve new model to get the elements up to speed
    new_model = Asap.Model(new_nodes, new_elements, new_loads)
    Asap.solve!(new_model, reprocess = true) # adjust all element lengths

    return new_model
end

"""
    update_load_values!(model::Asap.Model, params::SlabSizingParams, load_type::Symbol)

Updates the vertical load values (load.value[3]) for all loads in the model based on the specified load type from params.load_df.

# Arguments
- `model::Asap.Model`: The model containing loads to update
- `params::SlabSizingParams`: Parameters containing the load_df with load values
- `load_type::Symbol`: Which load value to use. Must be one of:
    - :factored_slab_load
    - :unfactored_slab_load  
    - :factored_applied_load
    - :unfactored_applied_load
    - :factored_beam_load
    - :unfactored_beam_load

"""
function update_load_values!(model::Asap.Model, params::SlabSizingParams; factored::Bool=true)
    loadid_index = _build_loadid_index(params)
    
    for load in model.loads
        # Point loads are mapped through load_df by loadID.
        if hasproperty(load, :loadID)
            new_load_value = get_new_load_value(params, load.loadID, factored=factored, loadid_index=loadid_index)
            if !isnothing(new_load_value)
                load.value = new_load_value
            end
        # Line loads (e.g., facade) do not carry loadID in this pipeline.
        elseif is_lineload(load)
            w_unfactored = !iszero(params.façade_load) ? abs(params.façade_load) :
                           (params.dead_factor > 0 ? abs(load.value[3]) / params.dead_factor : abs(load.value[3]))
            w_new = factored ? params.dead_factor * w_unfactored : w_unfactored
            load.value = [0, 0, -w_new]
        end
    end
end

function _build_loadid_index(params::SlabSizingParams)
    idx = Dict{Int, Int}()
    for (i, lid) in enumerate(params.load_df.loadID)
        idx[lid] = i
    end
    return idx
end

function get_new_load_value(params::SlabSizingParams, loadID::Int; factored::Bool=true, loadid_index::Union{Dict{Int,Int},Nothing}=nothing)
    slab_col = factored ? :factored_w_slab : :unfactored_w_slab
    applied_col = factored ? :factored_w_applied : :unfactored_w_applied 

    df_row = if !isnothing(loadid_index)
        get(loadid_index, loadID, nothing)
    else
        findfirst(==(loadID), params.load_df.loadID)
    end

    if !isnothing(df_row)
        slab_load = params.load_df[df_row, slab_col]
        applied_load = params.load_df[df_row, applied_col]
        total_load = slab_load + applied_load
        return -[0,0,total_load]
    else
        return nothing
    end
end

"""
    update_load_values_staged!(model, params; load_case)

Set model loads to a specific unfactored load case for staged deflection analysis.

# Load cases (all unfactored, for serviceability)
- `:slab_dead`  — slab self-weight only (bare steel stage, before concrete cures)
- `:sdl`        — superimposed dead load only (applied after composite action)
- `:live`       — live load only (applied after composite action)
- `:sdl_live`   — SDL + live (all post-composite loads)
- `:all`        — slab DL + SDL + live (total unfactored, equivalent to `factored=false`)

Façade line loads are treated as dead load (included in `:slab_dead` and `:all`).
"""
function update_load_values_staged!(model::Asap.Model, params::SlabSizingParams;
                                     load_case::Symbol=:all)
    loadid_index = _build_loadid_index(params)

    for load in model.loads
        if hasproperty(load, :loadID)
            df_row = get(loadid_index, load.loadID, nothing)
            isnothing(df_row) && continue

            w = if load_case == :slab_dead
                params.load_df[df_row, :unfactored_w_slab]
            elseif load_case == :sdl
                params.load_df[df_row, :unfactored_w_sdl]
            elseif load_case == :live
                params.load_df[df_row, :unfactored_w_live]
            elseif load_case == :sdl_live
                params.load_df[df_row, :unfactored_w_sdl] + params.load_df[df_row, :unfactored_w_live]
            else  # :all
                params.load_df[df_row, :unfactored_w_slab] +
                params.load_df[df_row, :unfactored_w_sdl] +
                params.load_df[df_row, :unfactored_w_live]
            end
            load.value = [0, 0, -w]

        elseif is_lineload(load)
            w_unfactored = !iszero(params.façade_load) ? abs(params.façade_load) :
                           (params.dead_factor > 0 ? abs(load.value[3]) / params.dead_factor : abs(load.value[3]))
            include_facade = load_case in (:slab_dead, :all)
            load.value = [0, 0, include_facade ? -w_unfactored : 0.0]
        end
    end
end

"""
    estimate_Ix(params::SlabSizingParams, beam_element::Element, beam_loads::Vector{Asap.AbstractLoad}, load_dictionary::Dict{String, Vector{Asap.AbstractLoad}}; material::Symbol=:steel_ksi)

Estimates required Ix using actual beam deflection from the current point-load pattern.
"""
function estimate_Ix(params::SlabSizingParams, beam_element::Element, beam_loads::Vector{Asap.AbstractLoad}, load_dictionary::Dict{String, Vector{Asap.AbstractLoad}}; material::Symbol=:steel_ksi)

    disp = ElementDisplacements(beam_element, beam_loads, resolution=200)
    δ = maximum(abs.(disp.ulocal[2, :])) / params.deflection_reduction_factor
    δ_max = beam_element.length / params.serviceability_lim
    Ix = beam_element.section.Ix

    return δ > 0 ? Ix * δ / δ_max : 0.0
end

"""
    estimate_δ(params::SlabSizingParams, beam_element::Element, beam_loads::Vector{Asap.AbstractLoad}, load_dictionary::Dict{String, Vector{Asap.AbstractLoad}}; material::Symbol=:steel_ksi)

Gets the maximum local vertical deflection from the current point-load pattern.
"""
function estimate_δ(params::SlabSizingParams, beam_element::Element, beam_loads::Vector{Asap.AbstractLoad}, load_dictionary::Dict{String, Vector{Asap.AbstractLoad}}; material::Symbol=:steel_ksi)

    disp = ElementDisplacements(beam_element, beam_loads, resolution=200)
    return maximum(abs.(disp.ulocal[2, :])) / params.deflection_reduction_factor
end

function sequential_search_sections(catalogue, objective, constraints, max_depth)
    valid_catalogue = filter(s -> s.d <= max_depth, catalogue)
    
    if isempty(valid_catalogue)
        throw(NoValidSectionsError("No sections found within max depth constraint"))
    end


    best_section = nothing
    best_objective = Inf

    # Iterate through all sections to find best valid one
    for section in valid_catalogue

        vars = get_geometry_vars(section)
        objective_val = objective(vars)

        # Check if section satisfies all constraints
        valid = true
        for constraint in constraints
            constraint_vector = constraint(vars)
            if any(x -> x > 0, constraint_vector)
                valid = false
                break
            end
        end

        # Update best section if valid and has better objective
        if valid && objective_val < best_objective
            best_section = section
            best_objective = objective_val
        end
    end

    if isnothing(best_section)
        throw(NoValidSectionsError("No valid sections found satisfying constraints"))
    end

    return best_objective, best_section
end

"""
    binary_search_sections(catalogue, objective, constraints, max_depth)

Performs binary search through sorted catalogue to find optimal section satisfying constraints.
Returns tuple of (best_objective_value, best_section).

# Arguments
- `catalogue`: Sorted collection of sections
- `objective`: Function that calculates objective value for a section
- `constraints`: Collection of constraint functions
- `max_depth`: Maximum allowable section depth

# Returns
- Tuple of (best_objective_value, best_section)
"""
function binary_search_sections(catalogue, objective, constraints, max_depth)
    # Filter out sections exceeding max depth and sort by objective
    valid_catalogue = filter(s -> s.d <= max_depth, catalogue)
    
    if isempty(valid_catalogue)
        throw(NoValidSectionsError("No sections found within max depth constraint"))
    end
    
    left = 1
    right = length(valid_catalogue)
    best_section = nothing
    best_objective = Inf
    
    while left <= right
        mid = Int(floor((left + right)/2))
        section = valid_catalogue[mid]
        vars = get_geometry_vars(section)
        objective_val = objective(vars)
        
        # Check constraints
        valid = true
        for constraint in constraints
            constraint_vector = constraint(vars)
            if any(x -> x > 0, constraint_vector)
                valid = false
                break
            end
        end
        
        if valid
            # Found valid section, look for better one in lower half
            if objective_val < best_objective
                best_section = section
                best_objective = objective_val
            end
            right = mid - 1
        else
            # Invalid section, look in upper half
            left = mid + 1
        end
    end
    
    if isnothing(best_section)
        throw(NoValidSectionsError("No valid sections found satisfying constraints"))
    end
    
    return best_objective, best_section
end

"""
    find_valid_sections(objective, constraints, max_depth)

Evaluates discrete sections for beam sizing; most naïve approch.
"""
function find_valid_sections(catalogue, objective, constraints, max_depth)
    valid_catalogue = filter(s -> s.d <= max_depth, catalogue)
    objective_vals, valid_sections = [], []

    for section in valid_catalogue
        if section.d > max_depth
            continue
        end

        vars = get_geometry_vars(section)
        objective_val = objective(vars)

        valid = true
        for constraint in constraints
            constraint_vector = constraint(vars)
            if any(x -> x > 0, constraint_vector)
                valid = false
                break
            end
        end

        if valid
            push!(objective_vals, objective_val)
            push!(valid_sections, section)
        end
    end

    return objective_vals, valid_sections
end