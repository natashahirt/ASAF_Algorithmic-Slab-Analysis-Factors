# ═══════════════════════════════════════════════════════════════════════════════
# Cached catalogue — built once from allW_imperial() at include-time.
# Both discrete paths (vector + MIP) and the continuous warm-start reuse this.
# ═══════════════════════════════════════════════════════════════════════════════

const _CACHED_CATALOGUE_DF = let
    catalogue = allW_imperial()
    n = length(catalogue)
    names  = Vector{String}(undef, n)
    h_v    = Vector{Float64}(undef, n)
    w_v    = Vector{Float64}(undef, n)
    tw_v   = Vector{Float64}(undef, n)
    tf_v   = Vector{Float64}(undef, n)
    A_v    = Vector{Float64}(undef, n)
    Ix_v   = Vector{Float64}(undef, n)
    Iy_v   = Vector{Float64}(undef, n)
    Iyc_v  = Vector{Float64}(undef, n)
    J_v    = Vector{Float64}(undef, n)
    Mn_v   = Vector{Float64}(undef, n)
    Vn_v   = Vector{Float64}(undef, n)

    for k in 1:n
        sec = I_symm(catalogue[k].d, catalogue[k].bf, catalogue[k].tw, catalogue[k].tf)
        names[k] = catalogue[k].name
        h_v[k]   = sec.h;   w_v[k]  = sec.w
        tw_v[k]  = sec.tw;  tf_v[k] = sec.tf
        A_v[k]   = sec.A;   Ix_v[k] = sec.Ix
        Iy_v[k]  = sec.Iy;  Iyc_v[k] = sec.I_yc
        J_v[k]   = sec.J
        Mn_v[k]  = sec.Mn;  Vn_v[k] = sec.Vn
    end

    DataFrame(name=names, h=h_v, w=w_v, tw=tw_v, tf=tf_v,
              A=A_v, Ix=Ix_v, Iy=Iy_v, I_yc=Iyc_v, J=J_v, Mn=Mn_v, Vn=Vn_v)
end

"""
    _get_catalogue_df(max_depth)

Return a depth-filtered view of the cached catalogue. Rows with h > max_depth
are excluded via a lightweight filter — no I_symm re-evaluation.
"""
function _get_catalogue_df(max_depth::Real)
    if iszero(max_depth) || isinf(max_depth)
        return _CACHED_CATALOGUE_DF
    end
    mask = _CACHED_CATALOGUE_DF.h .<= max_depth
    return _CACHED_CATALOGUE_DF[mask, :]
end

"""
    _active_beam_depth_limit(params) -> Float64

Return the active depth limit for beam sections. When assembly depth is in play,
`max_beam_depth` takes precedence over the raw `max_depth`; a zero or infinite
value means "unbounded".
"""
function _active_beam_depth_limit(params::SlabSizingParams)
    if !iszero(params.max_beam_depth) && !isinf(params.max_beam_depth)
        return Float64(params.max_beam_depth)
    end
    return Float64(params.max_depth)
end

"""
    _effective_composite_Ix(h, w, tw, tf, sp) -> Float64

Virtual-work effective composite Ix matching `get_I_composite_effective`,
but with all spatial constants (strip positions, effective widths, moment
weights) pre-computed in `sp`.  Used by both the NLP operator closure and
catalog fallback.

The effective Ix is:

    I_eff = sum_M2dx / Σ_k (M2dx_k / I_local_k)

where I_local_k = transformed-section composite Ix at strip k.
"""
function _effective_composite_Ix(h, w, tw, tf, sp::NamedTuple)
    t_slab = sp.t_slab
    n_mod  = sp.n_mod

    Ix_bare = Ix_I_asymm(h, w, w, tw, tf, tf)
    h_w = h - 2 * tf
    A_steel = 2 * w * tf + h_w * tw
    y_steel = h / 2
    y_conc  = h + t_slab / 2

    den = 0.0
    for k in 1:sp.n_strips
        b_tr_k     = sp.b_eff[k] / n_mod
        A_conc_tr_k = b_tr_k * t_slab
        A_total_k   = A_steel + A_conc_tr_k
        y_bar_k     = (A_steel * y_steel + A_conc_tr_k * y_conc) / A_total_k
        I_conc_own_k = b_tr_k * t_slab^3 / 12
        I_local_k = Ix_bare + A_steel * (y_steel - y_bar_k)^2 +
                    I_conc_own_k + A_conc_tr_k * (y_conc - y_bar_k)^2
        den += sp.M2dx[k] / I_local_k
    end

    den <= 0 && return Ix_bare
    return sp.sum_M2dx / den
end

"""
    _catalog_fallback_for_beam(;
        M_max, V_max, max_depth, Ix_bare_floor, Ix_comp_floor, comp_slab_params,
    ) -> Union{Vector{Float64}, Nothing}

Find the lightest W-shape in the cached catalogue that satisfies strength,
bare-Ix, and composite-Ix floor constraints for a single beam.  Returns
`[h, w, tw, tf]` or `nothing`.  The catalogue is pre-sorted by area
(ascending), so the first feasible hit is the lightest.
"""
function _catalog_fallback_for_beam(;
    M_max::Float64,
    V_max::Float64,
    max_depth::Float64,
    Ix_bare_floor::Float64=0.0,
    Ix_comp_floor::Float64=0.0,
    comp_slab_params::Union{Nothing, NamedTuple}=nothing,
)
    cat = _get_catalogue_df(max_depth)
    ϕ = 0.9
    for row in eachrow(cat)
        ϕ * row.Mn >= abs(M_max) || continue
        ϕ * row.Vn >= abs(V_max) || continue
        if Ix_bare_floor > 0
            row.Ix >= Ix_bare_floor || continue
        end
        if Ix_comp_floor > 0 && !isnothing(comp_slab_params)
            Ix_c = _effective_composite_Ix(row.h, row.w, row.tw, row.tf, comp_slab_params)
            Ix_c >= Ix_comp_floor || continue
        end
        return [row.h, row.w, row.tw, row.tf]
    end
    return nothing
end

"""
    _compute_staged_deflection_requirements!(
        params, beam_elements, beam_ids,
        bare_A, bare_Ix, sec_Iy, sec_J, comp_Ix;
        resolution=200,
    ) -> (Ix_live_req, Ix_coupled_comp, sw_coeff_240, δ_dead_total)

Shared staged-deflection helper for both MIP and NLP paths.  Performs two
global FE solves (slab DL on bare steel, SDL+LL on composite) and returns
per-beam Ix requirements:

- `Ix_live_req`    — composite Ix needed for L/360 live deflection
- `Ix_coupled_comp`— composite Ix needed for L/240 total (SDL+LL portion)
- `sw_coeff_240`   — self-weight coupling coefficient for the MIP constraint
- `δ_dead_total`   — dead-load deflection on bare steel (slab DL + beam SW)

Does **not** restore factored loads/sections — callers must do so.
"""
function _compute_staged_deflection_requirements!(
    params::SlabSizingParams,
    beam_elements::Vector{<:Asap.AbstractElement},
    beam_ids::Vector{Tuple{Int,Int}},
    bare_A::Vector{Float64},
    bare_Ix::Vector{Float64},
    sec_Iy::Vector{Float64},
    sec_J::Vector{Float64},
    comp_Ix::Vector{Float64};
    resolution::Int=200,
)
    n_beams = length(beam_elements)
    @assert length(beam_ids) == n_beams
    @assert length(bare_A) == n_beams
    @assert length(bare_Ix) == n_beams
    @assert length(sec_Iy) == n_beams
    @assert length(sec_J) == n_beams
    @assert length(comp_Ix) == n_beams

    E_s = steel_ksi.E
    ρ_steel = steel_ksi.ρ

    # ── Solve A: slab DL on bare steel ────────────────────────────────────
    for i in 1:n_beams
        beam_elements[i].section = Section(
            bare_A[i], E_s, steel_ksi.G, bare_Ix[i], sec_Iy[i], sec_J[i]
        )
    end
    update_load_values_staged!(params.model, params, load_case=:slab_dead)
    params.load_dictionary = get_load_dictionary_by_id(params.model)
    Asap.solve!(params.model, reprocess=true)

    δ_dead_bare = Vector{Float64}(undef, n_beams)
    for i in 1:n_beams
        disp = ElementDisplacements(
            beam_elements[i], params.load_dictionary[beam_ids[i]], resolution=resolution
        )
        δ_dead_bare[i] = maximum(abs.(disp.ulocal[2, :]))
    end

    # ── Solve B: SDL + LL on composite section ────────────────────────────
    for i in 1:n_beams
        beam_elements[i].section = Section(
            bare_A[i], E_s, steel_ksi.G, comp_Ix[i], sec_Iy[i], sec_J[i]
        )
    end
    update_load_values_staged!(params.model, params, load_case=:sdl_live)
    params.load_dictionary = get_load_dictionary_by_id(params.model)
    Asap.solve!(params.model, reprocess=true)

    δ_sdl_live_comp = Vector{Float64}(undef, n_beams)
    comp_Ix_used = Vector{Float64}(undef, n_beams)
    for i in 1:n_beams
        disp = ElementDisplacements(
            beam_elements[i], params.load_dictionary[beam_ids[i]], resolution=resolution
        )
        δ_sdl_live_comp[i] = maximum(abs.(disp.ulocal[2, :]))
        comp_Ix_used[i] = beam_elements[i].section.Ix
    end

    # ── Build per-beam requirements ───────────────────────────────────────
    Ix_live_req = Vector{Float64}(undef, n_beams)
    Ix_coupled_comp = Vector{Float64}(undef, n_beams)
    sw_coeff_240 = Vector{Float64}(undef, n_beams)
    δ_dead_total = Vector{Float64}(undef, n_beams)
    loadid_index_stg = _build_loadid_index(params)

    for i in 1:n_beams
        L_beam = beam_elements[i].length
        lim_360 = L_beam / 360.0
        lim_240 = L_beam / 240.0

        # Split SDL+LL deflection by load weight ratio
        w_sdl_sum = 0.0
        w_live_sum = 0.0
        for ld in params.load_dictionary[beam_ids[i]]
            if hasproperty(ld, :loadID)
                row = get(loadid_index_stg, ld.loadID, nothing)
                if !isnothing(row)
                    w_sdl_sum += params.load_df[row, :unfactored_w_sdl]
                    w_live_sum += params.load_df[row, :unfactored_w_live]
                end
            end
        end
        w_sum = w_sdl_sum + w_live_sum
        f_live_i = w_sum > 0 ? w_live_sum / w_sum : 0.5

        # Dead-load deflection on bare steel (slab DL + analytical beam SW)
        δ_sw_i = 5 * ρ_steel * bare_A[i] * L_beam^4 / (384 * E_s * bare_Ix[i])
        δ_dead_total[i] = δ_dead_bare[i] + δ_sw_i

        # L/360 live requirement (composite Ix)
        δ_live_i = δ_sdl_live_comp[i] * f_live_i
        Ix_live_req[i] = δ_live_i > 0 ? comp_Ix_used[i] * δ_live_i / lim_360 : 0.0

        # L/240 coupled requirement (composite Ix for SDL+LL portion)
        residual_240 = max(lim_240 - δ_dead_total[i], lim_240 * 0.1)
        Ix_coupled_comp[i] = δ_sdl_live_comp[i] > 0 ?
            comp_Ix_used[i] * δ_sdl_live_comp[i] / residual_240 : 0.0

        sw_coeff_240[i] = 5 * 240.0 * ρ_steel * L_beam^3 / (384 * E_s)
    end

    return Ix_live_req, Ix_coupled_comp, sw_coeff_240, δ_dead_total
end

"""
    process_discrete_beams_vector(params::SlabSizingParams)

Sizes beams using a vector sorting method. Not in parallel.
"""
function process_discrete_beams_vector(params::SlabSizingParams)

    model = params.model
    beam_elements = model.elements[:beam]
    best_sections = DataFrameRow[]
    active_max_depth = _active_beam_depth_limit(params)

    ϕ_b = 0.9
    ϕ_v = 0.9
    catalogue_df = _get_catalogue_df(active_max_depth)

    if isempty(catalogue_df)
        throw(NoValidSectionsError("No sections found within the active beam depth limit."))
    end

    for i in 1:lastindex(beam_elements)
        beam_id = get_element_id(beam_elements[i])
        beam_loads = params.load_dictionary[beam_id]
        beam_forces = InternalForces(beam_elements[i], beam_loads, resolution=200)
        M_max, V_max, x_max = maximum(abs.(beam_forces.My)), maximum(abs.(beam_forces.Vy)), maximum(abs.(beam_forces.x))

        feasible_sections = catalogue_df[(ϕ_b .* catalogue_df.Mn .>= M_max) .&
                                         (ϕ_v .* catalogue_df.Vn .>= V_max), :]

        if !isempty(feasible_sections)
            best_section = feasible_sections[argmin(feasible_sections.A), :]
            push!(best_sections, best_section)
        else
            throw(NoValidSectionsError("No valid sections found satisfying constraints"))
        end

        V = best_section.A * x_max

        push!(params.minimizers, [best_section.h, best_section.w, best_section.tw, best_section.tf])
        push!(params.minimums, V)
        push!(params.ids, best_section.name)
        push!(params.M_maxs, M_max)
        push!(params.V_maxs, V_max)
        push!(params.x_maxs, x_max)
    end

    println("Sized $(length(params.minimizers)) beams using vector sorting.")

    return params

end

"""
    integer_programming_discrete_beams(params::SlabSizingParams)

Sizes beams using integer programming with cached catalogue and threaded
composite Ix computation.
"""
function process_discrete_beams_integer(params::SlabSizingParams)

    n_max_sections = params.n_max_sections > 0 ? params.n_max_sections : nrow(_CACHED_CATALOGUE_DF)

    model = params.model
    beam_elements = model.elements[:beam]
    active_max_depth = _active_beam_depth_limit(params)

    catalogue_df = _get_catalogue_df(active_max_depth)
    n_cat = nrow(catalogue_df)

    if n_cat == 0
        throw(NoValidSectionsError("No sections found within the active beam depth limit."))
    end

    # Ensure factored beam self-weight is on the model so MIP demands include it
    refresh_factored_loads_with_beam_sw!(params)

    n_be = lastindex(beam_elements)
    beam_names = Vector{Tuple{Int,Int}}(undef, n_be)
    beam_M = Vector{Float64}(undef, n_be)
    beam_V = Vector{Float64}(undef, n_be)
    beam_x = Vector{Float64}(undef, n_be)

    for i in 1:n_be
        beam_id = get_element_id(beam_elements[i])
        beam_loads = params.load_dictionary[beam_id]
        beam_forces = InternalForces(beam_elements[i], beam_loads, resolution=200)
        beam_names[i] = beam_id
        beam_M[i] = maximum(abs.(beam_forces.My))
        beam_V[i] = maximum(abs.(beam_forces.Vy))
        beam_x[i] = maximum(abs.(beam_forces.x))
    end

    beams_df = DataFrame(name=beam_names, M_max=beam_M, V_max=beam_V, x_max=beam_x)

    if params.drawn
        element_ids = params.element_ids
        unique_ids = unique(element_ids)
        for element_id in unique_ids
            idxs = findall(x -> x == element_id, element_ids)
            max_M = maximum(beams_df.M_max[idxs])
            max_V = maximum(beams_df.V_max[idxs])
            max_x = maximum(beams_df.x_max[idxs])
            for i in idxs
                beams_df.M_max[i] = max_M
                beams_df.V_max[i] = max_V
                beams_df.x_max[i] = max_x
            end
        end
    end

    # ── Pre-filter by strength/depth BEFORE computing composite Ix ────────
    ϕ_b = 0.9
    ϕ_v = 0.9
    filtered_indices = Dict{Int, Vector{Int}}()
    for i in 1:n_be
        filtered_indices[i] = findall(
            (ϕ_b .* catalogue_df.Mn .>= beams_df.M_max[i]) .&
            (ϕ_v .* catalogue_df.Vn .>= beams_df.V_max[i])
        )
        if isempty(filtered_indices[i])
            throw(NoValidSectionsError(
                "No valid sections found for beam $i after strength/depth filtering " *
                "(Mu=$(round(beams_df.M_max[i], digits=2)), Vu=$(round(beams_df.V_max[i], digits=2)), " *
                "depth_limit=$(active_max_depth))."
            ))
        end
    end

    # ── Composite Ix via exact get_I_composite_effective (threaded) ───────
    E_s = steel_ksi.E
    loadid_index = _build_loadid_index(params)
    effective_Ix_arr = Vector{Vector{Float64}}(undef, n_be)

    if params.composite_action
        t_slab = params.slab_depth_in
        E_c = params.E_c

        beam_load_data = Vector{Tuple{Vector{Float64}, Vector{Float64}, Bool, Float64}}(undef, n_be)
        for i in 1:n_be
            element_loads_i = params.load_dictionary[get_element_id(beam_elements[i])]
            positions_i = Float64[]
            widths_i = Float64[]
            for ld in element_loads_i
                if hasproperty(ld, :loadID)
                    row = get(loadid_index, getproperty(ld, :loadID), nothing)
                    if !isnothing(row)
                        push!(positions_i, ld.position)
                        push!(widths_i, params.load_df[row, :trib_width])
                    end
                end
            end
            is_perim = i in params.i_perimeter
            L_beam = beam_elements[i].length
            beam_load_data[i] = (positions_i, widths_i, is_perim, L_beam)
        end

        Threads.@threads for i in 1:n_be
            positions_i, widths_i, is_perim, L_beam = beam_load_data[i]
            filt = filtered_indices[i]
            ix_vec = Vector{Float64}(undef, n_cat)
            ix_vec .= catalogue_df.Ix
            for j in filt
                ix_vec[j] = get_I_composite_effective(
                    catalogue_df.h[j], catalogue_df.w[j],
                    catalogue_df.tw[j], catalogue_df.tf[j],
                    t_slab, E_s, E_c, L_beam,
                    positions_i, widths_i;
                    is_perimeter=is_perim)
            end
            effective_Ix_arr[i] = ix_vec
        end
    else
        cat_ix = catalogue_df.Ix
        for i in 1:n_be
            effective_Ix_arr[i] = cat_ix
        end
    end

    effective_Ix = Dict{Int, Vector{Float64}}(i => effective_Ix_arr[i] for i in 1:n_be)

    ρ_steel = steel_ksi.ρ  # kip/in³

    has_staged = :unfactored_w_live in propertynames(params.load_df)
    use_staged_mip = params.deflection_limit && has_staged && params.composite_action

    # ── legacy (non-staged) deflection requirements ───────────────────────
    Ix_base_req = Dict{Int, Float64}()
    sw_coeff_req = Dict{Int, Float64}()

    if params.deflection_limit
        if !use_staged_mip
            update_load_values!(params.model, params, factored=false)
            params.load_dictionary = get_load_dictionary_by_id(params.model)
            Asap.solve!(params.model, reprocess=true)
            service_δ = Dict{Int, Float64}()
            for i in 1:lastindex(beam_elements)
                beam_id = get_element_id(beam_elements[i])
                element_loads_i = params.load_dictionary[beam_id]
                disp = ElementDisplacements(beam_elements[i], element_loads_i, resolution=200)
                service_δ[i] = maximum(abs.(disp.ulocal[2, :]))
            end

            update_load_values!(params.model, params, factored=true)
            params.load_dictionary = get_load_dictionary_by_id(params.model)
            Asap.solve!(params.model, reprocess=true)

            for i in 1:lastindex(beam_elements)
                L_beam = beam_elements[i].length
                δ_unfactored = get(service_δ, i, 0.0)
                Ix_initial = beam_elements[i].section.Ix
                δ_allowed = (L_beam / params.serviceability_lim) * params.deflection_reduction_factor
                Ix_base_req[i] = δ_unfactored > 0 ? Ix_initial * δ_unfactored / δ_allowed : 0.0
                sw_coeff_req[i] = 5 * params.serviceability_lim * ρ_steel * L_beam^3 / (384 * E_s * params.deflection_reduction_factor)
            end

            tol = 1e-9
            for i in 1:lastindex(beam_elements)
                feasible = [j for j in filtered_indices[i] if (effective_Ix[i][j] - sw_coeff_req[i] * catalogue_df.A[j] + tol) >= Ix_base_req[i]]
                if isempty(feasible)
                    throw(NoValidSectionsError("No valid sections found satisfying high-fidelity deflection constraint for beam $(i)."))
                end
                filtered_indices[i] = feasible
            end
        end
    end

    # ── apply outer-loop Ix floors from staged verification ─────────────
    if !isempty(params.min_Ix_comp) || !isempty(params.min_Ix_bare)
        tol = 1e-9
        for i in 1:lastindex(beam_elements)
            floor_comp = get(params.min_Ix_comp, i, 0.0)
            floor_bare = get(params.min_Ix_bare, i, 0.0)
            if floor_comp > 0 || floor_bare > 0
                feasible = [j for j in filtered_indices[i] if
                    effective_Ix[i][j] + tol >= floor_comp &&
                    catalogue_df.Ix[j] + tol >= floor_bare]
                if isempty(feasible)
                    throw(NoValidSectionsError(
                        "No section satisfies current Ix floors for beam $i " *
                        "(need comp=$(round(floor_comp, digits=1)), bare=$(round(floor_bare, digits=1)))."
                    ))
                end
                filtered_indices[i] = feasible
            end
            if params.deflection_limit && !use_staged_mip
                max_comp_avail = maximum(effective_Ix[i][j] for j in filtered_indices[i])
                Ix_base_req[i] = min(max(get(Ix_base_req, i, 0.0), floor_comp), max_comp_avail)
            end
        end
    end

    # Unify feasible catalog sets within collinear groups before the single MIP solve.
    groups = params.collinear_groups
    filtered_iter = deepcopy(filtered_indices)
    if !isempty(groups)
        for gid in unique(groups)
            members = findall(==(gid), groups)
            if length(members) > 1
                common = intersect([filtered_iter[m] for m in members]...)
                if isempty(common)
                    throw(NoValidSectionsError(
                        "Collinear group $gid has no common feasible sections under the current constraints."
                    ))
                end
                for m in members
                    filtered_iter[m] = common
                end
            end
        end
    end

    chosen_section_idx = Dict{Int, Int}()

    jump_model = JuMP.Model(Gurobi.Optimizer)
    JuMP.set_optimizer_attribute(jump_model, "OutputFlag", 1)
    JuMP.set_optimizer_attribute(jump_model, "MIPGap", 1e-4)
    JuMP.set_optimizer_attribute(jump_model, "NodeLimit", 10000)

    JuMP.@variable(jump_model, x[i=1:lastindex(beam_elements), j=filtered_iter[i]], Bin)
    JuMP.@constraint(jump_model, [i in 1:lastindex(beam_elements)],
        sum(x[i,j] for j in filtered_iter[i]) == 1)

    for i in 1:lastindex(beam_elements)
        JuMP.@constraint(jump_model, ϕ_b * sum(x[i,j] * catalogue_df.Mn[j] for j in filtered_iter[i]) >= beams_df.M_max[i])
        JuMP.@constraint(jump_model, ϕ_v * sum(x[i,j] * catalogue_df.Vn[j] for j in filtered_iter[i]) >= beams_df.V_max[i])
        if !iszero(active_max_depth) && !isinf(active_max_depth)
            JuMP.@constraint(jump_model, sum(x[i,j] * catalogue_df.h[j] for j in filtered_iter[i]) <= active_max_depth)
        end
        if params.deflection_limit && !use_staged_mip
            JuMP.@constraint(jump_model,
                sum(x[i,j] * (effective_Ix[i][j] - sw_coeff_req[i] * catalogue_df.A[j])
                    for j in filtered_iter[i]) >= Ix_base_req[i])
        end
    end

    if !isempty(groups)
        for gid in unique(groups)
            members = findall(==(gid), groups)
            if length(members) > 1
                leader = members[1]
                for k in 2:length(members)
                    follower = members[k]
                    for j in filtered_iter[leader]
                        JuMP.@constraint(jump_model, x[leader, j] == x[follower, j])
                    end
                end
            end
        end
    end

    JuMP.@variable(jump_model, y[j=1:lastindex(catalogue_df.name)], Bin)
    JuMP.@constraint(jump_model, [j=1:lastindex(catalogue_df.name)],
        sum(x[i,j] for i in 1:lastindex(beam_elements) if j in filtered_iter[i]) <= length(beam_elements) * y[j])
    JuMP.@constraint(jump_model, sum(y[j] for j in 1:lastindex(catalogue_df.name)) <= n_max_sections)

    JuMP.@objective(jump_model, Min,
        sum(sum(x[i,j] * catalogue_df.A[j] * beams_df.x_max[i] for j in filtered_iter[i])
            for i in 1:lastindex(beam_elements)))

    JuMP.optimize!(jump_model)

    term_status = JuMP.termination_status(jump_model)
    primal_status = JuMP.primal_status(jump_model)
    if !(term_status == JuMP.MOI.OPTIMAL || primal_status == JuMP.MOI.FEASIBLE_POINT)
        throw(NoValidSectionsError("Discrete MIP failed to produce a feasible solution (status=$(term_status))."))
    end
    term_status != JuMP.MOI.OPTIMAL && @warn "Discrete MIP returned a feasible but non-optimal solution." term_status primal_status

    for i in 1:lastindex(beam_elements)
        chosen_section_idx[i] = filtered_iter[i][argmax([JuMP.value(x[i,j]) for j in filtered_iter[i]])]
    end

    # Process results from final MIP solution
    for i in 1:lastindex(beam_elements)
        j = chosen_section_idx[i]
        best_section = catalogue_df[j, :]
        push!(params.minimizers, [best_section.h, best_section.w, best_section.tw, best_section.tf])
        push!(params.minimums, best_section.A * beams_df.x_max[i])
        push!(params.ids, best_section.name)
        push!(params.M_maxs, beams_df.M_max[i])
        push!(params.V_maxs, beams_df.V_max[i])
        push!(params.x_maxs, beams_df.x_max[i])
    end

    println("Sized $(length(params.minimizers)) beams using integer programming.")

    return params

end

"""
    _process_continuous_beams_parallel_legacy(...)

Legacy monolithic NLP beam sizer. Retained for A/B comparison; not called by
the default pipeline. Use `process_continuous_per_beam` instead.
"""
function _process_continuous_beams_parallel_legacy(params::SlabSizingParams;
                                                    initial_vars::Vector=[],
                                                    max_demand_iters::Int=3,
                                                    demand_tol::Float64=1e-3)

    model = params.model
    beam_elements = model.elements[:beam]
    n_beams = length(beam_elements)

    # ── initialisation ───────────────────────────────────────────────────
    min_h, min_w, min_tw, min_tf = params.minimum_continuous ? get_geometry_vars(W_imperial("W6X8.5")) : [0.01, 0.01, 0.001, 0.001]
    max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))

    if !iszero(params.max_beam_depth) && !isinf(params.max_beam_depth)
        max_h = min(params.max_beam_depth, max_h)
    end

    init_vars_array = Vector{Vector{Float64}}(undef, n_beams)
    for i in 1:n_beams
        if typeof(initial_vars) == Vector{String}
            init_vars_array[i] = get_geometry_vars(W_imperial(initial_vars[i]))
        elseif typeof(initial_vars) == Vector{Vector}
            init_vars_array[i] = initial_vars[i]
        elseif typeof(initial_vars) == Vector{Section}
            init_vars_array[i] = get_geometry_vars(initial_vars[i])
        else
            init_vars_array[i] = get_geometry_vars(W_imperial("W6X8.5"))
        end
    end

    ϕ_b = ϕ_v = 0.9
    Fy = steel_ksi.Fy
    E = steel_ksi.E
    loadid_index_cont = _build_loadid_index(params)

    # ── NL helper functions (defined once, registered per JuMP model) ────
    @inline _constraint_area(h, w, tw, tf) = -A_I_asymm(h, w, w, tw, tf, tf)
    @inline _constraint_height(h, tf) = 2 * tf - h
    @inline function _constraint_web_flange_ratio(h, w, tw, tf)
        h_w = h - 2 * tf
        return tw * h_w / (tf * w) - 10
    end
    @inline _constraint_web_slenderness(h, tw, tf) = (h - 2 * tf) / tw - 260
    @inline function _constraint_flange_ratio_min(h, w, tw, tf)
        I_yc = parallel_axis_I(tf, w, 0)
        Iy = Iy_I_asymm(h, w, w, tw, tf, tf)
        return 0.1 - I_yc/Iy
    end
    @inline function _constraint_flange_ratio_max(h, w, tw, tf)
        I_yc = parallel_axis_I(tf, w, 0)
        Iy = Iy_I_asymm(h, w, w, tw, tf, tf)
        return I_yc/Iy - 0.9
    end
    @inline function _Vn_I_symm(h, w, tw, tf, Fy, E)
        h_w = h - 2 * tf
        Aw = h * tw
        kv = 5.34
        lim_rolled = 2.24 * sqrt(E / Fy)
        lim_inelastic = 1.10 * sqrt(kv * E / Fy)
        ratio = h_w / tw
        Cv1 = ratio <= lim_rolled ? 1.0 :
              (ratio <= lim_inelastic ? 1.0 : lim_inelastic / ratio)
        return 0.6 * Fy * Aw * Cv1
    end

    @inline function _Mn_LTB_I_symm(h, w, tw, tf, Fy, E, Lb)
        A  = A_I_asymm(h, w, w, tw, tf, tf)
        Ix = Ix_I_asymm(h, w, w, tw, tf, tf)
        Iy = Iy_I_asymm(h, w, w, tw, tf, tf)
        J  = J_I_symm(h, w, tw, tf)
        Zx = Zx_I_asymm(h, w, w, tw, tf, tf)
        Sx = Ix / (h / 2)
        ry = sqrt(Iy / max(A, 1e-10))
        ho = h - tf
        Cw = Iy * ho^2 / 4
        rts = Sx > 0 ? sqrt(sqrt(Iy * Cw) / Sx) : 1e-10

        Mp = Fy * Zx
        Mn = Mp

        if Lb > 0
            c  = 1.0
            Cb = 1.0
            Lp = 1.76 * ry * sqrt(E / Fy)

            jc_term = (J * c) / max(Sx * ho, 1e-10)
            Lr = 1.95 * rts * (E / (0.7 * Fy)) *
                 sqrt(jc_term + sqrt(jc_term^2 + 6.76 * (0.7 * Fy / E)^2))

            if Lb > Lr
                lb_rts = Lb / max(rts, 1e-10)
                F_cr = Cb * π^2 * E / lb_rts^2 *
                       sqrt(1 + 0.078 * jc_term * lb_rts^2)
                Mn = min(Mn, F_cr * Sx)
            elseif Lb > Lp
                Mn_LTB = Cb * (Mp - (Mp - 0.7 * Fy * Sx) * ((Lb - Lp) / max(Lr - Lp, 1e-10)))
                Mn = min(Mn, min(Mn_LTB, Mp))
            end
        end
        return Mn
    end

    # ── precompute composite load data (geometry only, once) ─────────────
    local beam_load_data_cont
    if params.composite_action
        beam_load_data_cont = Vector{Tuple{Vector{Float64}, Vector{Float64}, Bool, Float64}}(undef, n_beams)
        for i in 1:n_beams
            element_loads_i = params.load_dictionary[get_element_id(beam_elements[i])]
            positions_i = Float64[]
            widths_i = Float64[]
            for ld in element_loads_i
                if hasproperty(ld, :loadID)
                    row = get(loadid_index_cont, getproperty(ld, :loadID), nothing)
                    if !isnothing(row)
                        push!(positions_i, ld.position)
                        push!(widths_i, params.load_df[row, :trib_width])
                    end
                end
            end
            is_perim = i in params.i_perimeter
            L_beam = beam_elements[i].length
            beam_load_data_cont[i] = (positions_i, widths_i, is_perim, L_beam)
        end
    end

    prev_obj = Inf
    current_vars = copy(init_vars_array)

    # ── Fix A: set FE model sections to warm-start values before demand loop ──
    if !isempty(initial_vars)
        for i in 1:n_beams
            hi, wi, twi, tfi = current_vars[i]
            A  = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
            Ix = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
            Iy = Iy_I_asymm(hi, wi, wi, twi, tfi, tfi)
            J  = J_I_symm(hi, wi, twi, tfi)
            beam_elements[i].section = Section(A, steel_ksi.E, steel_ksi.G, Ix, Iy, J)
        end
        Asap.solve!(params.model, reprocess=true)
        params.load_dictionary = get_load_dictionary_by_id(params.model)
    end

    # ── Fix C: warm-start areas for per-iteration upper-bound constraint ──
    warm_start_areas = Vector{Float64}(undef, n_beams)
    has_warm_start = !isempty(initial_vars)
    if has_warm_start
        for i in 1:n_beams
            hi, wi, twi, tfi = current_vars[i]
            warm_start_areas[i] = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
        end
    end

    # ── Fix B: track best solution across demand iterations ──
    best_vars = copy(current_vars)
    best_obj  = Inf

    for demand_iter in 1:max_demand_iters
        println("─── demand iteration $demand_iter / $max_demand_iters ───")

        # Factored slab/applied loads + beam self-weight line loads for strength demands
        update_load_values!(params.model, params, factored=true)
        sync_beam_selfweight_lineloads!(params, factored=true, vars=current_vars)
        Asap.solve!(params.model, reprocess=true)
        params.load_dictionary = get_load_dictionary_by_id(params.model)

        # ── compute demands from current stiffnesses ─────────────────────
        beam_ids = Vector{Tuple{Int,Int}}(undef, n_beams)
        M_maxs = Vector{Float64}(undef, n_beams)
        V_maxs = Vector{Float64}(undef, n_beams)
        x_maxs = Vector{Float64}(undef, n_beams)

        for i in 1:n_beams
            beam_ids[i] = get_element_id(beam_elements[i])
            beam_loads = params.load_dictionary[beam_ids[i]]
            M_maxs[i], V_maxs[i], x_maxs[i] = internalforces_M_V(beam_elements[i], beam_loads, resolution=200)
        end

        beams_df = DataFrame(
            name = beam_ids,
            M_max = M_maxs,
            V_max = V_maxs,
            x_max = x_maxs
        )

        # ── recompute deflection requirements from current sections ──────
        # Uses staged construction logic when load decomposition is available:
        #   Solve A: slab DL + facade on bare steel  →  δ_dead_bare
        #   Beam SW: analytical 5wL⁴/(384EI) on bare steel  →  δ_sw
        #   Solve B: SDL + LL on composite steel  →  δ_sdl_live (split by load ratio)
        # Then enforce L/360 live and L/240 total.
        # Falls back to single-solve L/serviceability_lim when staged columns are absent.
        local Ix_req_live_vec, Ix_coupled_comp_vec, Ix_req_bare_vec, sw_coeff_vec
        local use_staged_cont
        if params.deflection_limit
            has_staged_cont = :unfactored_w_live in propertynames(params.load_df)
            use_staged_cont = has_staged_cont && params.composite_action

            if use_staged_cont
                E_s = steel_ksi.E

                # Compute fresh section properties from current_vars
                fresh_A  = Vector{Float64}(undef, n_beams)
                fresh_Ix = Vector{Float64}(undef, n_beams)
                fresh_Iy = Vector{Float64}(undef, n_beams)
                fresh_J  = Vector{Float64}(undef, n_beams)
                comp_Ix = Vector{Float64}(undef, n_beams)
                for i in 1:n_beams
                    hi, wi, twi, tfi = current_vars[i]
                    fresh_A[i]  = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
                    fresh_Ix[i] = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
                    fresh_Iy[i] = Iy_I_asymm(hi, wi, wi, twi, tfi, tfi)
                    fresh_J[i]  = J_I_symm(hi, wi, twi, tfi)
                    if params.composite_action && params.slab_depth_in > 0
                        positions_i, widths_i, is_perim, L_beam = beam_load_data_cont[i]
                        if !isempty(positions_i)
                            comp_Ix[i] = get_I_composite_effective(hi, wi, twi, tfi,
                                params.slab_depth_in, E_s, params.E_c, L_beam,
                                positions_i, widths_i; is_perimeter=is_perim)
                        else
                            comp_Ix[i] = fresh_Ix[i]
                        end
                    else
                        comp_Ix[i] = fresh_Ix[i]
                    end
                end

                Ix_req_live_vec, Ix_coupled_comp_vec, _, _ =
                    _compute_staged_deflection_requirements!(
                        params, beam_elements, beam_ids,
                        fresh_A, fresh_Ix, fresh_Iy, fresh_J, comp_Ix;
                        resolution=200,
                    )

                # Restore bare-steel sections and factored loads (including beam self-weight)
                for i in 1:n_beams
                    beam_elements[i].section = Section(fresh_A[i], E_s, steel_ksi.G,
                        fresh_Ix[i], fresh_Iy[i], fresh_J[i])
                end
                refresh_factored_loads_with_beam_sw!(params)

                Ix_req_bare_vec = zeros(n_beams)
                sw_coeff_vec = zeros(n_beams)
            else
                # Legacy single-solve fallback (non-staged)
                update_load_values!(params.model, params, factored=false)
                params.load_dictionary = get_load_dictionary_by_id(params.model)
                Asap.solve!(params.model, reprocess=true)

                service_δ = Dict{Int, Float64}()
                for i in 1:n_beams
                    beam_id = get_element_id(beam_elements[i])
                    element_loads_i = params.load_dictionary[beam_id]
                    disp = ElementDisplacements(beam_elements[i], element_loads_i, resolution=200)
                    service_δ[i] = maximum(abs.(disp.ulocal[2, :]))
                end

                update_load_values!(params.model, params, factored=true)
                params.load_dictionary = get_load_dictionary_by_id(params.model)
                Asap.solve!(params.model, reprocess=true)

                ρ_steel_cont = steel_ksi.ρ
                E_s = steel_ksi.E
                Ix_req_live_vec       = zeros(n_beams)
                Ix_coupled_comp_vec = zeros(n_beams)
                Ix_req_bare_vec       = Vector{Float64}(undef, n_beams)
                sw_coeff_vec          = Vector{Float64}(undef, n_beams)
                for i in 1:n_beams
                    L_beam = beam_elements[i].length
                    δ_unfactored = get(service_δ, i, 0.0)
                    Ix_initial = beam_elements[i].section.Ix
                    δ_allowed = (L_beam / params.serviceability_lim) * params.deflection_reduction_factor
                    Ix_req_bare_vec[i] = δ_unfactored > 0 ? Ix_initial * δ_unfactored / δ_allowed : 0.0
                    sw_coeff_vec[i] = 5 * params.serviceability_lim * ρ_steel_cont * L_beam^3 / (384 * E_s * params.deflection_reduction_factor)
                end
            end
        end

        # ── compute composite ratio from current vars (exact per beam) ───
        local composite_ratio_vec
        needs_composite_ratio = params.composite_action &&
            (params.deflection_limit || !isempty(params.min_Ix_comp))
        if needs_composite_ratio
            t_slab = params.slab_depth_in
            E_c = params.E_c
            composite_ratio_vec = ones(n_beams)
            for i in 1:n_beams
                hi, wi, twi, tfi = current_vars[i]
                Ix_bare = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
                positions_i, widths_i, is_perim, L_beam = beam_load_data_cont[i]
                if !isempty(positions_i) && Ix_bare > 0
                    Ix_comp = get_I_composite_effective(hi, wi, twi, tfi,
                        t_slab, E, E_c, L_beam, positions_i, widths_i;
                        is_perimeter=is_perim)
                    composite_ratio_vec[i] = Ix_comp / Ix_bare
                end
            end
        end

        # ── build JuMP model ─────────────────────────────────────────────
        nlopt_algorithms = Dict(
            :MMA    => :LD_MMA,
            :SLSQP  => :LD_SLSQP,
            :CCSAQ  => :LD_CCSAQ,
            :COBYLA => :LN_COBYLA,
        )

        if haskey(nlopt_algorithms, params.nlp_solver)
            jump_model = JuMP.Model(NLopt.Optimizer)
            JuMP.set_attribute(jump_model, "algorithm", nlopt_algorithms[params.nlp_solver])
            JuMP.set_attribute(jump_model, "xtol_rel", 1e-4)
            JuMP.set_attribute(jump_model, "ftol_rel", 1e-4)
            JuMP.set_attribute(jump_model, "maxeval", 5000)
            JuMP.set_attribute(jump_model, "constrtol_abs", 1e-6)
        elseif params.nlp_solver == :Ipopt
            jump_model = JuMP.Model(Ipopt.Optimizer)
            JuMP.set_optimizer_attribute(jump_model, "print_level", demand_iter == 1 ? 1 : 0)
            JuMP.set_optimizer_attribute(jump_model, "warm_start_init_point", "yes")
            JuMP.set_optimizer_attribute(jump_model, "mu_init", 1e-4)
            JuMP.set_optimizer_attribute(jump_model, "tol", 1e-6)
            JuMP.set_optimizer_attribute(jump_model, "max_iter", 5000)
            JuMP.set_optimizer_attribute(jump_model, "bound_push", 1e-8)
            JuMP.set_optimizer_attribute(jump_model, "bound_frac", 1e-8)
            JuMP.set_optimizer_attribute(jump_model, "mu_strategy", "adaptive")
        else
            error("Unknown NLP solver: $(params.nlp_solver). " *
                  "Supported: :MMA, :SLSQP, :CCSAQ, :COBYLA, :Ipopt")
        end

        JuMP.@variable(jump_model, min_h <= h[i=1:n_beams] <= max_h, start=current_vars[i][1])
        JuMP.@variable(jump_model, min_w <= w[i=1:n_beams] <= max_w, start=current_vars[i][2])
        JuMP.@variable(jump_model, min_tw <= tw[i=1:n_beams] <= max_tw, start=current_vars[i][3])
        JuMP.@variable(jump_model, min_tf <= tf[i=1:n_beams] <= max_tf, start=current_vars[i][4])

        JuMP.@operator(jump_model, op_A_I_asymm, 6, A_I_asymm)
        JuMP.@operator(jump_model, op_Iy_I_asymm, 6, Iy_I_asymm)
        JuMP.@operator(jump_model, op_parallel_axis_I, 3, parallel_axis_I)
        JuMP.@operator(jump_model, op_Zx_I_asymm, 6, Zx_I_asymm)
        JuMP.@operator(jump_model, op_Ix_I_asymm, 6, Ix_I_asymm)
        JuMP.@operator(jump_model, op_constraint_area, 4, _constraint_area)
        JuMP.@operator(jump_model, op_constraint_height, 2, _constraint_height)
        JuMP.@operator(jump_model, op_constraint_web_flange_ratio, 4, _constraint_web_flange_ratio)
        JuMP.@operator(jump_model, op_constraint_web_slenderness, 3, _constraint_web_slenderness)
        JuMP.@operator(jump_model, op_constraint_flange_ratio_min, 4, _constraint_flange_ratio_min)
        JuMP.@operator(jump_model, op_constraint_flange_ratio_max, 4, _constraint_flange_ratio_max)
        JuMP.@operator(jump_model, op_Vn_I_symm, 6, _Vn_I_symm)
        JuMP.@operator(jump_model, op_Mn_LTB_I_symm, 7, _Mn_LTB_I_symm)

        # Unbraced length per beam: composite → 0 (slab braces top flange),
        # non-composite → full beam span.
        Lb_vec = [params.composite_action ? 0.0 : beam_elements[i].length for i in 1:n_beams]

        JuMP.@objective(jump_model, Min, sum(
            op_A_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) * beams_df.x_max[i] for i in 1:n_beams)
        )

        # Fix C: prevent optimizer from exceeding warm-start objective
        # Recompute bound each iteration using current beams_df spans so the
        # constraint stays consistent with the JuMP objective after demand
        # redistribution.
        if has_warm_start
            ws_obj = sum(warm_start_areas[i] * beams_df.x_max[i] for i in 1:n_beams)
            JuMP.@constraint(jump_model, sum(
                op_A_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) * beams_df.x_max[i]
                for i in 1:n_beams) <= ws_obj * 1.01)
        end

        for i in 1:n_beams
            Lb_i = Lb_vec[i]
            JuMP.@constraint(jump_model, op_constraint_area(h[i], w[i], tw[i], tf[i]) <= 0)
            JuMP.@constraint(jump_model, op_constraint_height(h[i], tf[i]) <= 0)
            JuMP.@constraint(jump_model, op_constraint_web_flange_ratio(h[i], w[i], tw[i], tf[i]) <= 0)
            JuMP.@constraint(jump_model, op_constraint_web_slenderness(h[i], tw[i], tf[i]) <= 0)
            JuMP.@constraint(jump_model, w[i] / (2 * tf[i]) <= 0.38 * sqrt(E/Fy))
            JuMP.@constraint(jump_model, (h[i] - 2 * tf[i]) / tw[i] <= 3.76 * sqrt(E/Fy))
            JuMP.@constraint(jump_model,
                abs(beams_df.M_max[i]) - ϕ_b * op_Mn_LTB_I_symm(h[i], w[i], tw[i], tf[i], Fy, E, Lb_i) <= 0)
            JuMP.@constraint(jump_model,
                abs(beams_df.V_max[i]) - ϕ_v * op_Vn_I_symm(h[i], w[i], tw[i], tf[i], Fy, E) <= 0)
        end

        # ── deflection constraints ───────────────────────────────────────
        # Staged: L/360 live (on composite Ix), L/240 total-comp (on composite Ix),
        #         L/240 dead (on bare Ix minus self-weight deduction).
        # Legacy: single Ix_req on effective Ix (composite or bare).
        if params.deflection_limit
            if use_staged_cont
                for i in 1:n_beams
                    if needs_composite_ratio && composite_ratio_vec[i] > 1.0
                        ratio_i = composite_ratio_vec[i]
                        JuMP.@constraint(jump_model,
                            Ix_req_live_vec[i] - ratio_i * op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
                        JuMP.@constraint(jump_model,
                            Ix_coupled_comp_vec[i] - ratio_i * op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
                    else
                        JuMP.@constraint(jump_model,
                            Ix_req_live_vec[i] - op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
                        JuMP.@constraint(jump_model,
                            Ix_coupled_comp_vec[i] - op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
                    end
                end
            else
                for i in 1:n_beams
                    sw_c = sw_coeff_vec[i]
                    if needs_composite_ratio && composite_ratio_vec[i] > 1.0
                        ratio_i = composite_ratio_vec[i]
                        JuMP.@constraint(jump_model,
                            Ix_req_bare_vec[i] - (ratio_i * op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) - sw_c * op_A_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i])) <= 0)
                    else
                        JuMP.@constraint(jump_model,
                            Ix_req_bare_vec[i] - (op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) - sw_c * op_A_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i])) <= 0)
                    end
                end
            end
        end

        # ── collinearity constraints: same section for beams in each group ─
        groups = params.collinear_groups
        if !isempty(groups)
            for gid in unique(groups)
                members = findall(==(gid), groups)
                if length(members) > 1
                    leader = members[1]
                    for k in 2:length(members)
                        fol = members[k]
                        JuMP.@constraint(jump_model, h[fol] == h[leader])
                        JuMP.@constraint(jump_model, w[fol] == w[leader])
                        JuMP.@constraint(jump_model, tw[fol] == tw[leader])
                        JuMP.@constraint(jump_model, tf[fol] == tf[leader])
                    end
                end
            end
        end

        # ── outer-loop Ix floor constraints from staged verification ──────
        for i in 1:n_beams
            floor_comp = get(params.min_Ix_comp, i, 0.0)
            floor_bare = get(params.min_Ix_bare, i, 0.0)

            if floor_comp > 0 && needs_composite_ratio && composite_ratio_vec[i] > 1.0
                ratio_i = composite_ratio_vec[i]
                JuMP.@constraint(jump_model,
                    floor_comp - ratio_i * op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
            elseif floor_comp > 0
                JuMP.@constraint(jump_model,
                    floor_comp - op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
            end

            if floor_bare > 0
                JuMP.@constraint(jump_model,
                    floor_bare - op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
            end
        end

        # ── solve ────────────────────────────────────────────────────────
        JuMP.optimize!(jump_model)

        status = JuMP.termination_status(jump_model)
        if status ∉ (JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.OPTIMAL,
                      JuMP.MOI.ALMOST_LOCALLY_SOLVED, JuMP.MOI.ALMOST_OPTIMAL)
            @warn "Continuous optimization did not converge " *
                  "(status: $status, demand iter $demand_iter/$max_demand_iters). " *
                  "Using best solution found so far."
            break
        end

        optimal_h = JuMP.value.(h)
        optimal_w = JuMP.value.(w)
        optimal_tw = JuMP.value.(tw)
        optimal_tf = JuMP.value.(tf)

        current_vars = [[optimal_h[i], optimal_w[i], optimal_tw[i], optimal_tf[i]] for i in 1:n_beams]
        area = [A_I_symm(current_vars[i]...) for i in 1:n_beams]
        cur_obj = sum(area[i] * beams_df.x_max[i] for i in 1:n_beams)

        println("  objective = $cur_obj")
        rel_change = abs(cur_obj - prev_obj) / max(abs(prev_obj), 1e-12)
        println("  relative change = $rel_change")

        # Fix B: keep the best solution seen so far
        if cur_obj < best_obj
            best_obj  = cur_obj
            best_vars = deepcopy(current_vars)
        end

        if demand_iter > 1 && rel_change < demand_tol
            println("  demand loop converged after $demand_iter iterations.")
            break
        end
        prev_obj = cur_obj

        # ── update FE model sections (skip full I_symm — inline props) ───
        if demand_iter < max_demand_iters
            for i in 1:n_beams
                hi, wi, twi, tfi = current_vars[i]
                A  = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
                Ix = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
                Iy = Iy_I_asymm(hi, wi, wi, twi, tfi, tfi)
                J  = J_I_symm(hi, wi, twi, tfi)
                beam_elements[i].section = Section(A, steel_ksi.E, steel_ksi.G, Ix, Iy, J)
            end
            Asap.solve!(params.model, reprocess=true)
            params.load_dictionary = get_load_dictionary_by_id(params.model)
        end
    end

    # Fix B: use the best solution found across all demand iterations
    current_vars = best_vars
    # Update FE model to reflect best solution for downstream postprocessing
    for i in 1:n_beams
        hi, wi, twi, tfi = current_vars[i]
        A  = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
        Ix = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
        Iy = Iy_I_asymm(hi, wi, wi, twi, tfi, tfi)
        J  = J_I_symm(hi, wi, twi, tfi)
        beam_elements[i].section = Section(A, steel_ksi.E, steel_ksi.G, Ix, Iy, J)
    end
    update_load_values!(params.model, params, factored=true)
    sync_beam_selfweight_lineloads!(params, factored=true, vars=current_vars)
    Asap.solve!(params.model, reprocess=true)
    params.load_dictionary = get_load_dictionary_by_id(params.model)

    # ── store results ────────────────────────────────────────────────────
    minimizer = current_vars
    area = [A_I_symm(minimizer[i]...) for i in 1:n_beams]

    final_M = Vector{Float64}(undef, n_beams)
    final_V = Vector{Float64}(undef, n_beams)
    final_x = Vector{Float64}(undef, n_beams)
    for i in 1:n_beams
        beam_id = get_element_id(beam_elements[i])
        beam_loads = params.load_dictionary[beam_id]
        final_M[i], final_V[i], final_x[i] = internalforces_M_V(beam_elements[i], beam_loads, resolution=200)
    end

    # Strict final gate: reject continuous candidate if final staged/strength checks fail.
    final_feasible = true
    tol_feas = 1e-6
    strength_fail_limit = 1.02

    has_staged_final = params.deflection_limit &&
        params.composite_action &&
        (:unfactored_w_live in propertynames(params.load_df))

    if has_staged_final
        E_s = steel_ksi.E
        fresh_A = Vector{Float64}(undef, n_beams)
        fresh_Ix = Vector{Float64}(undef, n_beams)
        fresh_Iy = Vector{Float64}(undef, n_beams)
        fresh_J = Vector{Float64}(undef, n_beams)
        comp_Ix = Vector{Float64}(undef, n_beams)

        for i in 1:n_beams
            hi, wi, twi, tfi = current_vars[i]
            fresh_A[i] = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
            fresh_Ix[i] = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
            fresh_Iy[i] = Iy_I_asymm(hi, wi, wi, twi, tfi, tfi)
            fresh_J[i] = J_I_symm(hi, wi, twi, tfi)

            if params.slab_depth_in > 0
                positions_i, widths_i, is_perim, L_beam = beam_load_data_cont[i]
                comp_Ix[i] = isempty(positions_i) ? fresh_Ix[i] : get_I_composite_effective(
                    hi, wi, twi, tfi, params.slab_depth_in, E_s, params.E_c, L_beam,
                    positions_i, widths_i; is_perimeter=is_perim
                )
            else
                comp_Ix[i] = fresh_Ix[i]
            end
        end

        final_beam_ids = [get_element_id(be) for be in beam_elements]
        Ix_live_final, Ix_coupled_final, _, _ = _compute_staged_deflection_requirements!(
            params, beam_elements, final_beam_ids,
            fresh_A, fresh_Ix, fresh_Iy, fresh_J, comp_Ix;
            resolution=200,
        )

        # Restore bare sections and factored loads after staged requirement recompute
        for i in 1:n_beams
            beam_elements[i].section = Section(
                fresh_A[i], E_s, steel_ksi.G, fresh_Ix[i], fresh_Iy[i], fresh_J[i]
            )
        end
        refresh_factored_loads_with_beam_sw!(params)

        for i in 1:n_beams
            req_comp = max(Ix_live_final[i], Ix_coupled_final[i], get(params.min_Ix_comp, i, 0.0))
            req_bare = get(params.min_Ix_bare, i, 0.0)
            if comp_Ix[i] + tol_feas < req_comp || fresh_Ix[i] + tol_feas < req_bare
                final_feasible = false
                break
            end
        end
    elseif params.deflection_limit
        for i in 1:n_beams
            hi, wi, twi, tfi = current_vars[i]
            A_i = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
            Ix_i = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
            if (Ix_i - sw_coeff_vec[i] * A_i) + tol_feas < Ix_req_bare_vec[i]
                final_feasible = false
                break
            end
        end
    end

    if final_feasible
        for i in 1:n_beams
            sec = I_symm(current_vars[i]...)
            Mn = get_ϕMn(sec)
            Vn = get_ϕVn(sec)
            if (Mn > 0 && final_M[i] / Mn > strength_fail_limit) ||
               (Vn > 0 && final_V[i] / Vn > strength_fail_limit)
                final_feasible = false
                break
            end
        end
    end

    if !final_feasible
        @warn "Rejecting continuous solution: final staged/strength feasibility check failed."
        params.minimizers = Vector{Float64}[]
        params.minimums = Float64[]
        params.ids = String[]
        params.M_maxs = Float64[]
        params.V_maxs = Float64[]
        params.x_maxs = Float64[]
        return params
    end

    minimum = [area[i] * final_x[i] for i in 1:n_beams]

    println("Sized $(length(minimizer)) beams using simultaneous optimization with demand iteration.")

    params.minimizers = minimizer
    params.minimums = minimum
    params.ids = string.(round.(area, digits=2))
    params.M_maxs = final_M
    params.V_maxs = final_V
    params.x_maxs = final_x

    return params

end

# ═══════════════════════════════════════════════════════════════════════════════
# Per-beam NLP refinement using frozen MIP demands
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _solve_beam_nlp(;
        M_max, V_max, x_max, Lb, init_vars,
        min_Ix_bare_floor, min_Ix_comp_floor, comp_slab_params,
        bounds, nlp_solver,
    ) -> Union{Vector{Float64}, Nothing}

Solve a 4-variable NLP (h, w, tw, tf) for a single beam with frozen demands
from the MIP.  Deflection is constrained via two channels:

- **Bare-Ix floor** — dead-load deflection on bare steel (L/240 budget)
- **Composite-Ix floor** — live L/360 and coupled L/240 via the full spatially-varying
  effective composite Ix (virtual-work weighted over strip positions), computed from
  the NLP's own geometry variables

Returns `[h, w, tw, tf]` or `nothing` on failure.
"""
function _solve_beam_nlp(;
    M_max::Float64,
    V_max::Float64,
    x_max::Float64,
    Lb::Float64,
    init_vars::Vector{Float64},
    min_Ix_bare_floor::Float64=0.0,
    min_Ix_comp_floor::Float64=0.0,
    comp_slab_params::Union{Nothing, NamedTuple}=nothing,
    bounds::NamedTuple,
    nlp_solver::Symbol=:MMA,
)
    ϕ_b = ϕ_v = 0.9
    Fy = steel_ksi.Fy
    E  = steel_ksi.E

    # NL constraint helpers (closures over nothing — pure geometry)
    _ca(h, w, tw, tf) = -A_I_asymm(h, w, w, tw, tf, tf)
    _ch(h, tf) = 2 * tf - h
    function _cwfr(h, w, tw, tf)
        h_w = h - 2 * tf
        return tw * h_w / (tf * w) - 10
    end
    _cws(h, tw, tf) = (h - 2 * tf) / tw - 260

    function _Vn(h, w, tw, tf, Fy_v, E_v)
        h_w = h - 2 * tf
        Aw = h * tw
        kv = 5.34
        lim_rolled = 2.24 * sqrt(E_v / Fy_v)
        lim_inelastic = 1.10 * sqrt(kv * E_v / Fy_v)
        ratio = h_w / tw
        Cv1 = ratio <= lim_rolled ? 1.0 :
              (ratio <= lim_inelastic ? 1.0 : lim_inelastic / ratio)
        return 0.6 * Fy_v * Aw * Cv1
    end

    function _Mn(h, w, tw, tf, Fy_v, E_v, Lb_v)
        A  = A_I_asymm(h, w, w, tw, tf, tf)
        Ix = Ix_I_asymm(h, w, w, tw, tf, tf)
        Iy = Iy_I_asymm(h, w, w, tw, tf, tf)
        J  = J_I_symm(h, w, tw, tf)
        Zx = Zx_I_asymm(h, w, w, tw, tf, tf)
        Sx = Ix / (h / 2)
        ry = sqrt(Iy / max(A, 1e-10))
        ho = h - tf
        Cw = Iy * ho^2 / 4
        rts = Sx > 0 ? sqrt(sqrt(Iy * Cw) / Sx) : 1e-10
        Mp = Fy_v * Zx
        Mn_val = Mp
        if Lb_v > 0
            c  = 1.0
            Cb = 1.0
            Lp = 1.76 * ry * sqrt(E_v / Fy_v)
            jc_term = (J * c) / max(Sx * ho, 1e-10)
            Lr = 1.95 * rts * (E_v / (0.7 * Fy_v)) *
                 sqrt(jc_term + sqrt(jc_term^2 + 6.76 * (0.7 * Fy_v / E_v)^2))
            if Lb_v > Lr
                lb_rts = Lb_v / max(rts, 1e-10)
                F_cr = Cb * π^2 * E_v / lb_rts^2 *
                       sqrt(1 + 0.078 * jc_term * lb_rts^2)
                Mn_val = min(Mn_val, F_cr * Sx)
            elseif Lb_v > Lp
                Mn_LTB = Cb * (Mp - (Mp - 0.7 * Fy_v * Sx) * ((Lb_v - Lp) / max(Lr - Lp, 1e-10)))
                Mn_val = min(Mn_val, min(Mn_LTB, Mp))
            end
        end
        return Mn_val
    end

    # Build JuMP model
    nlopt_algorithms = Dict(
        :MMA    => :LD_MMA,
        :SLSQP  => :LD_SLSQP,
        :CCSAQ  => :LD_CCSAQ,
        :COBYLA => :LN_COBYLA,
    )

    if haskey(nlopt_algorithms, nlp_solver)
        jm = JuMP.Model(NLopt.Optimizer)
        JuMP.set_attribute(jm, "algorithm", nlopt_algorithms[nlp_solver])
        JuMP.set_attribute(jm, "xtol_rel", 1e-6)
        JuMP.set_attribute(jm, "ftol_rel", 1e-6)
        JuMP.set_attribute(jm, "maxeval", 3000)
        JuMP.set_attribute(jm, "constrtol_abs", 1e-6)
    elseif nlp_solver == :Ipopt
        jm = JuMP.Model(Ipopt.Optimizer)
        JuMP.set_optimizer_attribute(jm, "print_level", 0)
        JuMP.set_optimizer_attribute(jm, "warm_start_init_point", "yes")
        JuMP.set_optimizer_attribute(jm, "mu_init", 1e-4)
        JuMP.set_optimizer_attribute(jm, "tol", 1e-6)
        JuMP.set_optimizer_attribute(jm, "max_iter", 3000)
        JuMP.set_optimizer_attribute(jm, "bound_push", 1e-8)
        JuMP.set_optimizer_attribute(jm, "bound_frac", 1e-8)
        JuMP.set_optimizer_attribute(jm, "mu_strategy", "adaptive")
    else
        error("Unknown NLP solver: $nlp_solver")
    end

    JuMP.@variable(jm, bounds.min_h <= h <= bounds.max_h, start=init_vars[1])
    JuMP.@variable(jm, bounds.min_w <= w <= bounds.max_w, start=init_vars[2])
    JuMP.@variable(jm, bounds.min_tw <= tw <= bounds.max_tw, start=init_vars[3])
    JuMP.@variable(jm, bounds.min_tf <= tf <= bounds.max_tf, start=init_vars[4])

    # Register nonlinear operators on this model
    JuMP.@operator(jm, op_A, 6, A_I_asymm)
    JuMP.@operator(jm, op_Ix, 6, Ix_I_asymm)
    JuMP.@operator(jm, op_ca, 4, _ca)
    JuMP.@operator(jm, op_ch, 2, _ch)
    JuMP.@operator(jm, op_cwfr, 4, _cwfr)
    JuMP.@operator(jm, op_cws, 3, _cws)
    JuMP.@operator(jm, op_Vn, 6, _Vn)
    JuMP.@operator(jm, op_Mn, 7, _Mn)

    # Objective: minimize cross-sectional area
    JuMP.@objective(jm, Min, op_A(h, w, w, tw, tf, tf))

    # AISC geometry constraints
    JuMP.@constraint(jm, op_ca(h, w, tw, tf) <= 0)
    JuMP.@constraint(jm, op_ch(h, tf) <= 0)
    JuMP.@constraint(jm, op_cwfr(h, w, tw, tf) <= 0)
    JuMP.@constraint(jm, op_cws(h, tw, tf) <= 0)
    JuMP.@constraint(jm, w / (2 * tf) <= 0.38 * sqrt(E / Fy))
    JuMP.@constraint(jm, (h - 2 * tf) / tw <= 3.76 * sqrt(E / Fy))

    # Strength constraints
    JuMP.@constraint(jm, abs(M_max) - ϕ_b * op_Mn(h, w, tw, tf, Fy, E, Lb) <= 0)
    JuMP.@constraint(jm, abs(V_max) - ϕ_v * op_Vn(h, w, tw, tf, Fy, E) <= 0)

    # Deflection: bare-Ix lower bound (dead-load on bare steel)
    if min_Ix_bare_floor > 0
        JuMP.@constraint(jm, min_Ix_bare_floor - op_Ix(h, w, w, tw, tf, tf) <= 0)
    end

    # Deflection: composite-Ix lower bound (live L/360, total L/240 on composite section).
    # Uses the full spatially-varying effective composite Ix (virtual-work weighted
    # over strip positions) — same formula as get_I_composite_effective.
    if min_Ix_comp_floor > 0 && !isnothing(comp_slab_params)
        _sp = comp_slab_params

        _Ix_comp_fn = (h_v, w_v, tw_v, tf_v) -> begin
            _effective_composite_Ix(h_v, w_v, tw_v, tf_v, _sp)
        end

        JuMP.@operator(jm, op_Ix_comp, 4, _Ix_comp_fn)
        JuMP.@constraint(jm, min_Ix_comp_floor - op_Ix_comp(h, w, tw, tf) <= 0)
    end

    JuMP.optimize!(jm)

    status = JuMP.termination_status(jm)
    if status ∉ (JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.OPTIMAL,
                  JuMP.MOI.ALMOST_LOCALLY_SOLVED, JuMP.MOI.ALMOST_OPTIMAL)
        return nothing
    end

    return [JuMP.value(h), JuMP.value(w), JuMP.value(tw), JuMP.value(tf)]
end

"""
    process_continuous_per_beam(params; mip_result)

Per-beam NLP refinement using frozen MIP demands.

For each beam (or collinear group leader), solves a 4-variable NLP minimizing
cross-sectional area subject to:

- **Strength**: MIP-frozen M_max/V_max demands with AISC geometry constraints
- **Deflection** (two channels):
    - *Bare-Ix floor* — dead-load deflection on bare steel ≤ L/240 budget
    - *Composite-Ix floor* — L/360 live and L/240 total requirements, constrained
      directly on the NLP's own spatially-varying effective composite Ix (same
      virtual-work formula as `get_I_composite_effective`)
    - Legacy non-staged deflection (when staged loads are absent)

The outer staged-deflection ratchet loop in `optimal_beamsizer` is skipped for
the NLP path.  Falls back to catalog discrete section or MIP section when the
NLP is infeasible.
"""
function process_continuous_per_beam(params::SlabSizingParams;
                                     mip_result::Union{SlabSizingParams, Nothing}=nothing)
    if isnothing(mip_result) || isempty(mip_result.minimizers)
        @warn "process_continuous_per_beam: no MIP result available; returning empty."
        params.minimizers = Vector{Float64}[]
        params.minimums = Float64[]
        params.ids = String[]
        params.M_maxs = Float64[]
        params.V_maxs = Float64[]
        params.x_maxs = Float64[]
        return params
    end

    model = params.model
    beam_elements = model.elements[:beam]
    n_beams = length(beam_elements)

    @assert length(mip_result.minimizers) == n_beams "MIP result beam count mismatch"
    @assert length(mip_result.M_maxs) == n_beams

    # Frozen demands from MIP
    mip_M = copy(mip_result.M_maxs)
    mip_V = copy(mip_result.V_maxs)
    mip_x = copy(mip_result.x_maxs)
    mip_vars = deepcopy(mip_result.minimizers)

    # Staged Ix requirements and composite ratios are always computed from MIP
    # geometry (single-pass — no outer-loop iteration for NLP).
    staging_vars = mip_vars

    E_s = steel_ksi.E

    # Variable bounds (same as legacy NLP)
    min_h, min_w, min_tw, min_tf = params.minimum_continuous ?
        get_geometry_vars(W_imperial("W6X8.5")) : [0.01, 0.01, 0.001, 0.001]
    max_h, max_w, max_tw, max_tf = get_geometry_vars(W_imperial("W43X335"))
    if !iszero(params.max_beam_depth) && !isinf(params.max_beam_depth)
        max_h = min(params.max_beam_depth, max_h)
    end
    var_bounds = (min_h=min_h, max_h=max_h, min_w=min_w, max_w=max_w,
                  min_tw=min_tw, max_tw=max_tw, min_tf=min_tf, max_tf=max_tf)

    # Unbraced length per beam
    Lb_vec = [params.composite_action ? 0.0 : beam_elements[i].length for i in 1:n_beams]

    # Compute composite load data for composite Ix ratios
    loadid_index = _build_loadid_index(params)
    local beam_load_data
    if params.composite_action
        beam_load_data = Vector{Tuple{Vector{Float64}, Vector{Float64}, Bool, Float64}}(undef, n_beams)
        for i in 1:n_beams
            element_loads_i = params.load_dictionary[get_element_id(beam_elements[i])]
            positions_i = Float64[]
            widths_i = Float64[]
            for ld in element_loads_i
                if hasproperty(ld, :loadID)
                    row = get(loadid_index, getproperty(ld, :loadID), nothing)
                    if !isnothing(row)
                        push!(positions_i, ld.position)
                        push!(widths_i, params.load_df[row, :trib_width])
                    end
                end
            end
            is_perim = i in params.i_perimeter
            L_beam = beam_elements[i].length
            beam_load_data[i] = (positions_i, widths_i, is_perim, L_beam)
        end
    end

    # Pre-compute spatially-varying composite slab constants per beam for the
    # NLP constraint.  This mirrors get_I_composite_effective: aggregated strip
    # positions, per-strip effective widths, and virtual-work moment weights
    # M̄(x)²·dx are frozen so the NLP closure only evaluates the steel-geometry-
    # dependent terms at each strip.
    comp_slab_params_vec = Vector{Union{Nothing, NamedTuple}}(undef, n_beams)
    fill!(comp_slab_params_vec, nothing)
    if params.composite_action && params.slab_depth_in > 0
        n_mod = E_s / params.E_c
        for i in 1:n_beams
            positions_i, widths_i, is_perim, L_beam = beam_load_data[i]
            if !isempty(widths_i)
                agg_pos, agg_w = _aggregate_strips(positions_i, widths_i)
                n_pts = length(agg_pos)

                x_abs = agg_pos .* L_beam
                M_bar = x_abs .* (L_beam .- x_abs)

                perm = sortperm(x_abs)
                x_sorted = x_abs[perm]
                dx = zeros(n_pts)
                if n_pts == 1
                    dx[perm[1]] = L_beam
                else
                    for (idx, si) in enumerate(perm)
                        if idx == 1
                            dx[si] = (x_sorted[2] - x_sorted[1]) / 2 + x_sorted[1]
                        elseif idx == n_pts
                            dx[si] = (L_beam - x_sorted[end]) + (x_sorted[end] - x_sorted[end-1]) / 2
                        else
                            dx[si] = (x_sorted[idx+1] - x_sorted[idx-1]) / 2
                        end
                    end
                end

                M2dx     = [M_bar[k]^2 * dx[k] for k in 1:n_pts]
                b_eff_k  = [get_b_eff(L_beam, agg_w[k]; is_perimeter=is_perim) for k in 1:n_pts]
                sum_M2dx = sum(M2dx)

                if sum_M2dx > 0
                    comp_slab_params_vec[i] = (
                        t_slab   = params.slab_depth_in,
                        n_mod    = n_mod,
                        n_strips = n_pts,
                        M2dx     = M2dx,
                        b_eff    = b_eff_k,
                        sum_M2dx = sum_M2dx,
                    )
                end
            end
        end
    end

    # ── Compute Ix floors per beam from MIP staging geometry ──────────
    #
    # The postprocessor checks staged deflection (L/360 live, L/240 total) via
    # FE solves.  Requirements are split into two channels:
    #
    #   Bare-Ix floor  — dead-load deflection on bare steel (L/240 budget)
    #   Comp-Ix floor  — live L/360 and coupled L/240 (composite Ix, constrained directly)
    #   Legacy         — non-staged bare-Ix requirement (when staged loads absent)

    use_staged = params.deflection_limit &&
        params.composite_action &&
        (:unfactored_w_live in propertynames(params.load_df))

    # Section properties from MIP geometry for staged FE solves
    bare_A  = Vector{Float64}(undef, n_beams)
    bare_Ix = Vector{Float64}(undef, n_beams)
    sec_Iy  = Vector{Float64}(undef, n_beams)
    sec_J   = Vector{Float64}(undef, n_beams)
    comp_Ix = Vector{Float64}(undef, n_beams)

    for i in 1:n_beams
        hi, wi, twi, tfi = staging_vars[i]
        bare_A[i]  = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
        bare_Ix[i] = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
        sec_Iy[i]  = Iy_I_asymm(hi, wi, wi, twi, tfi, tfi)
        sec_J[i]   = J_I_symm(hi, wi, twi, tfi)
        if params.composite_action && params.slab_depth_in > 0
            positions_i, widths_i, is_perim, L_beam = beam_load_data[i]
            if !isempty(positions_i)
                comp_Ix[i] = get_I_composite_effective(hi, wi, twi, tfi,
                    params.slab_depth_in, E_s, params.E_c, L_beam,
                    positions_i, widths_i; is_perimeter=is_perim)
            else
                comp_Ix[i] = bare_Ix[i]
            end
        else
            comp_Ix[i] = bare_Ix[i]
        end
    end

    # Run staged FE solves and build requirement vectors
    env_Ix_bare_floor = zeros(n_beams)
    env_Ix_comp_floor = zeros(n_beams)

    if params.deflection_limit && use_staged
        beam_ids = [get_element_id(be) for be in beam_elements]
        Ix_req_live_vec, Ix_coupled_comp_vec, _, δ_dead_total_staging =
            _compute_staged_deflection_requirements!(
                params, beam_elements, beam_ids,
                bare_A, bare_Ix, sec_Iy, sec_J, comp_Ix;
                resolution=200,
            )

        for i in 1:n_beams
            beam_elements[i].section = Section(
                bare_A[i], E_s, steel_ksi.G, bare_Ix[i], sec_Iy[i], sec_J[i])
        end
        refresh_factored_loads_with_beam_sw!(params)

        for i in 1:n_beams
            # Dead-load deflection on bare steel → bare-Ix floor
            L_beam = beam_elements[i].length
            lim_240 = L_beam / 240.0
            if lim_240 > 0 && bare_Ix[i] > 0 && δ_dead_total_staging[i] > 0
                env_Ix_bare_floor[i] = bare_Ix[i] * δ_dead_total_staging[i] / lim_240
            end

            # Live L/360 and coupled L/240 → composite-Ix floor (direct, no ratio conversion)
            comp_floor = 0.0
            if Ix_req_live_vec[i] > 0
                comp_floor = max(comp_floor, Ix_req_live_vec[i])
            end
            if Ix_coupled_comp_vec[i] > 0
                comp_floor = max(comp_floor, Ix_coupled_comp_vec[i])
            end

            if isnothing(comp_slab_params_vec[i]) && comp_floor > 0
                # No slab load data for this beam (composite ratio = 1.0);
                # the composite requirement reduces to a bare-Ix requirement.
                env_Ix_bare_floor[i] = max(env_Ix_bare_floor[i], comp_floor)
            else
                env_Ix_comp_floor[i] = comp_floor
            end
        end
    elseif params.deflection_limit
        # Legacy non-staged: single unfactored deflection solve
        for i in 1:n_beams
            beam_elements[i].section = Section(
                bare_A[i], E_s, steel_ksi.G, bare_Ix[i], sec_Iy[i], sec_J[i])
        end
        update_load_values!(params.model, params, factored=false)
        params.load_dictionary = get_load_dictionary_by_id(params.model)
        Asap.solve!(params.model, reprocess=true)

        for i in 1:n_beams
            beam_id = get_element_id(beam_elements[i])
            disp = ElementDisplacements(beam_elements[i],
                params.load_dictionary[beam_id], resolution=200)
            δ_unfactored = maximum(abs.(disp.ulocal[2, :]))
            Ix_initial = beam_elements[i].section.Ix
            L_beam = beam_elements[i].length
            δ_allowed = (L_beam / params.serviceability_lim) * params.deflection_reduction_factor
            if δ_unfactored > 0
                env_Ix_bare_floor[i] = Ix_initial * δ_unfactored / δ_allowed
            end
        end

        update_load_values!(params.model, params, factored=true)
        params.load_dictionary = get_load_dictionary_by_id(params.model)
        Asap.solve!(params.model, reprocess=true)
    end

    # ── Collinear grouping: envelope demands across group members ─────────
    groups = params.collinear_groups
    beam_to_leader = collect(1:n_beams)
    group_members = Dict{Int, Vector{Int}}()
    if !isempty(groups)
        for gid in unique(groups)
            members = findall(==(gid), groups)
            if length(members) > 1
                leader = members[1]
                for m in members
                    beam_to_leader[m] = leader
                end
                group_members[leader] = members
            end
        end
    end

    env_M = copy(mip_M)
    env_V = copy(mip_V)
    env_x = copy(mip_x)

    for (leader, members) in group_members
        env_M[leader] = maximum(mip_M[m] for m in members)
        env_V[leader] = maximum(mip_V[m] for m in members)
        env_x[leader] = maximum(mip_x[m] for m in members)
        env_Ix_bare_floor[leader] = maximum(env_Ix_bare_floor[m] for m in members)
        env_Ix_comp_floor[leader] = maximum(env_Ix_comp_floor[m] for m in members)
        # For composite slab params, use the member whose strip layout produces the
        # lowest effective composite Ix at the MIP geometry (most demanding).
        worst_Ix = Inf
        worst_idx = leader
        h_l, w_l, tw_l, tf_l = staging_vars[leader]
        for m in members
            sp = comp_slab_params_vec[m]
            if !isnothing(sp)
                Ix_c = _effective_composite_Ix(h_l, w_l, tw_l, tf_l, sp)
                if Ix_c < worst_Ix
                    worst_Ix = Ix_c
                    worst_idx = m
                end
            end
        end
        comp_slab_params_vec[leader] = comp_slab_params_vec[worst_idx]
    end

    # Solve per beam (or per group leader).  Deflection requirements are split
    # into a bare-Ix floor (dead load) and a composite-Ix floor (live/total),
    # both constrained directly in the NLP.  When the NLP is infeasible, a
    # catalog lookup finds the lightest discrete W-shape satisfying the floors,
    # then warm-starts a retry.
    result_vars = Vector{Vector{Float64}}(undef, n_beams)
    n_fallback_mip = 0
    n_fallback_catalog = 0
    n_catalog_retry_ok = 0

    leaders_to_solve = unique(beam_to_leader)
    max_depth_for_catalog = var_bounds.max_h

    for leader in leaders_to_solve
        sol = _solve_beam_nlp(
            M_max  = env_M[leader],
            V_max  = env_V[leader],
            x_max  = env_x[leader],
            Lb     = Lb_vec[leader],
            init_vars = staging_vars[leader],
            min_Ix_bare_floor  = env_Ix_bare_floor[leader],
            min_Ix_comp_floor  = env_Ix_comp_floor[leader],
            comp_slab_params   = comp_slab_params_vec[leader],
            bounds     = var_bounds,
            nlp_solver = params.nlp_solver,
        )

        # Catalog-lookup fallback: find the lightest discrete section that
        # satisfies the current floors, then retry NLP with it as warm start.
        if isnothing(sol)
            cat_vars = _catalog_fallback_for_beam(
                M_max  = env_M[leader],
                V_max  = env_V[leader],
                max_depth = max_depth_for_catalog,
                Ix_bare_floor  = env_Ix_bare_floor[leader],
                Ix_comp_floor  = env_Ix_comp_floor[leader],
                comp_slab_params = comp_slab_params_vec[leader],
            )
            if !isnothing(cat_vars)
                retry_sol = _solve_beam_nlp(
                    M_max  = env_M[leader],
                    V_max  = env_V[leader],
                    x_max  = env_x[leader],
                    Lb     = Lb_vec[leader],
                    init_vars = cat_vars,
                    min_Ix_bare_floor  = env_Ix_bare_floor[leader],
                    min_Ix_comp_floor  = env_Ix_comp_floor[leader],
                    comp_slab_params   = comp_slab_params_vec[leader],
                    bounds     = var_bounds,
                    nlp_solver = params.nlp_solver,
                )
                if !isnothing(retry_sol)
                    sol = retry_sol
                    n_catalog_retry_ok += 1
                else
                    sol = cat_vars
                    n_fallback_catalog += 1
                end
            end
        end

        if isnothing(sol)
            sol = mip_vars[leader]
            n_fallback_mip += 1
        end

        members = get(group_members, leader, [leader])
        for m in members
            result_vars[m] = sol
        end
    end

    n_leaders = length(leaders_to_solve)
    if n_catalog_retry_ok > 0
        println("  Per-beam NLP: $n_catalog_retry_ok / $n_leaders beam(s) recovered via catalog warm-start.")
    end
    if n_fallback_catalog > 0
        println("  Per-beam NLP: $n_fallback_catalog / $n_leaders beam(s) using catalog discrete section (NLP retry failed).")
    end
    if n_fallback_mip > 0
        println("  Per-beam NLP: $n_fallback_mip / $n_leaders beam(s) fell back to MIP section (no catalog match).")
    end

    function _apply_result_sections!(vars_set::Vector{Vector{Float64}})
        for i in 1:n_beams
            hi, wi, twi, tfi = vars_set[i]
            A  = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
            Ix = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
            Iy = Iy_I_asymm(hi, wi, wi, twi, tfi, tfi)
            J  = J_I_symm(hi, wi, twi, tfi)
            beam_elements[i].section = Section(A, E_s, steel_ksi.G, Ix, Iy, J)
        end
    end

    function _verify_final_candidate!(vars_set::Vector{Vector{Float64}})
        params.minimizers = deepcopy(vars_set)
        finalize_beam_selfweight_factored_demands!(params)

        final_M = copy(params.M_maxs)
        final_V = copy(params.V_maxs)
        final_x = copy(params.x_maxs)
        failing_beams = Set{Int}()
        failing_leaders = Set{Int}()
        staged_n_viol = 0

        for i in 1:n_beams
            sec = I_symm(vars_set[i]...)
            ϕMn = get_ϕMn(sec)
            ϕVn = get_ϕVn(sec)
            if (ϕMn > 0 && final_M[i] / ϕMn > 1.02) ||
               (ϕVn > 0 && final_V[i] / ϕVn > 1.02)
                push!(failing_beams, i)
                push!(failing_leaders, beam_to_leader[i])
            end
        end

        if params.deflection_limit && use_staged
            staged_n_viol, new_Ix_comp, new_Ix_bare = _verify_staged_deflection(params)
            for i in union(Set(keys(new_Ix_comp)), Set(keys(new_Ix_bare)))
                push!(failing_beams, i)
                push!(failing_leaders, beam_to_leader[i])
            end
        elseif params.deflection_limit
            _apply_result_sections!(vars_set)
            update_load_values!(params.model, params, factored=false)
            sync_beam_selfweight_lineloads!(params, factored=false, vars=vars_set)
            params.load_dictionary = get_load_dictionary_by_id(params.model)
            Asap.solve!(params.model, reprocess=true)

            for i in 1:n_beams
                beam_id = get_element_id(beam_elements[i])
                disp = ElementDisplacements(
                    beam_elements[i], params.load_dictionary[beam_id], resolution=200
                )
                δ_unfactored = maximum(abs.(disp.ulocal[2, :]))
                δ_allowed = (beam_elements[i].length / params.serviceability_lim) *
                            params.deflection_reduction_factor
                if δ_unfactored > δ_allowed + 1e-6
                    push!(failing_beams, i)
                    push!(failing_leaders, beam_to_leader[i])
                end
            end
            refresh_factored_loads_with_beam_sw!(params)
        end

        return failing_beams, failing_leaders, final_M, final_V, final_x, staged_n_viol
    end

    final_M = copy(mip_M)
    final_V = copy(mip_V)
    final_x = copy(mip_x)
    final_reverted_leaders = Set{Int}()
    final_staged_viol = 0

    for verify_iter in 1:3
        failing_beams, failing_leaders, M_chk, V_chk, x_chk, staged_chk =
            _verify_final_candidate!(result_vars)
        final_M, final_V, final_x = M_chk, V_chk, x_chk
        final_staged_viol = staged_chk

        if isempty(failing_beams)
            break
        end

        union!(final_reverted_leaders, failing_leaders)
        for leader in failing_leaders
            members = get(group_members, leader, [leader])
            for m in members
                result_vars[m] = mip_vars[m]
            end
        end
    end

    if !isempty(final_reverted_leaders)
        println("  Per-beam NLP: reverted $(length(final_reverted_leaders)) leader/group(s) to MIP after final verification.")
        _, _, final_M, final_V, final_x, final_staged_viol = _verify_final_candidate!(result_vars)
    end

    # Leave the global model in the accepted final factored state.
    _apply_result_sections!(result_vars)
    update_load_values!(params.model, params, factored=true)
    sync_beam_selfweight_lineloads!(params, factored=true, vars=result_vars)
    Asap.solve!(params.model, reprocess=true)
    params.load_dictionary = get_load_dictionary_by_id(params.model)

    area = [A_I_symm(result_vars[i]...) for i in 1:n_beams]

    params.minimizers = result_vars
    params.minimums   = [area[i] * final_x[i] for i in 1:n_beams]
    params.ids        = string.(round.(area, digits=2))
    params.M_maxs     = final_M
    params.V_maxs     = final_V
    params.x_maxs     = final_x
    params.staged_n_violations = final_staged_viol
    params.staged_converged = final_staged_viol == 0

    total_obj = sum(area[i] * final_x[i] for i in 1:n_beams)
    mip_obj = sum(A_I_symm(mip_vars[i]...) * mip_x[i] for i in 1:n_beams)
    pct = mip_obj > 0 ? round(100 * (1 - total_obj / mip_obj), digits=2) : 0.0
    println("Sized $n_beams beams via per-beam NLP refinement ($(pct)% lighter than MIP).")

    return params
end

function internalforces_M_V(beam_element::Element, beam_loads::Vector{Asap.AbstractLoad}; resolution::Int=200)

    uglobal = [beam_element.nodeStart.displacement; beam_element.nodeEnd.displacement]
    release = AsapToolkit.get_release_type(beam_element)
    r2dof = AsapToolkit.release2DOF[release]
    Flocal = (beam_element.R * beam_element.K * uglobal) .* r2dof
    Vystart =  Flocal[2]
    Mystart =  Flocal[6]
    xinc = range(0, beam_element.length; length=resolution)
    R = beam_element.R[1:3, 1:3]
    L = beam_element.length

    moment_function = AsapToolkit.MPointLoad[release]
    shear_function = AsapToolkit.VPointLoad[release]

    My = Vector{Float64}(undef, lastindex(xinc))
    Vy = Vector{Float64}(undef, lastindex(xinc))

    My .= Vystart .* xinc .- Mystart
    Vy .= Vystart

    for load in beam_loads
        px, py, pz = (R * load.value) .* [1, -1, -1]
        if is_pointload(load)
            My .+= moment_function.(py, L, xinc, load.position)
            Vy .+= shear_function.(py, L, xinc, load.position)
        elseif is_lineload(load)
            My .+= AsapToolkit.MLine.(Ref(beam_element), py, L, xinc)
            Vy .+= AsapToolkit.VLine.(Ref(beam_element), py, L, xinc)
        end
    end

    return maximum(abs.(My)), maximum(abs.(Vy)), L

end

function internalforces_M_V_all(beam_element::Element, beam_loads::Vector{Asap.AbstractLoad}; resolution::Int=200)

    uglobal = [beam_element.nodeStart.displacement; beam_element.nodeEnd.displacement]
    release = AsapToolkit.get_release_type(beam_element)
    r2dof = AsapToolkit.release2DOF[release]
    Flocal = (beam_element.R * beam_element.K * uglobal) .* r2dof
    Vystart =  Flocal[2]
    Mystart =  Flocal[6]
    xinc = range(0, beam_element.length; length=resolution)
    R = beam_element.R[1:3, 1:3]
    L = beam_element.length

    moment_function = AsapToolkit.MPointLoad[release]
    shear_function = AsapToolkit.VPointLoad[release]

    My = Vector{Float64}(undef, lastindex(xinc))
    Vy = Vector{Float64}(undef, lastindex(xinc))

    My .= Vystart .* xinc .- Mystart
    Vy .= Vystart

    for load in beam_loads
        px, py, pz = (R * load.value) .* [1, -1, -1]
        if is_pointload(load)
            My .+= moment_function.(py, L, xinc, load.position)
            Vy .+= shear_function.(py, L, xinc, load.position)
        elseif is_lineload(load)
            My .+= AsapToolkit.MLine.(Ref(beam_element), py, L, xinc)
            Vy .+= AsapToolkit.VLine.(Ref(beam_element), py, L, xinc)
        end
    end

    return abs.(My), abs.(Vy), collect(xinc)

end