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
    process_discrete_beams_vector(params::SlabSizingParams)

Sizes beams using a vector sorting method. Not in parallel.
"""
function process_discrete_beams_vector(params::SlabSizingParams)

    model = params.model
    beam_elements = model.elements[:beam]
    best_sections = DataFrameRow[]

    ϕ_b = 0.9
    ϕ_v = 0.9
    catalogue_df = _get_catalogue_df(params.max_depth)

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

    catalogue_df = _get_catalogue_df(params.max_depth)
    n_cat = nrow(catalogue_df)

    # Build beams DataFrame column-wise
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

    # ── staged deflection requirement computation ─────────────────────────
    # Extracted so the iterative re-solve loop can call it with updated sections.

    function _compute_staged_Ix_requirements(section_indices::Dict{Int, Int})
        Ix_live_req       = Dict{Int, Float64}()
        Ix_coupled_comp   = Dict{Int, Float64}()
        Ix_coupled_bare   = Dict{Int, Float64}()
        sw_coeff_240      = Dict{Int, Float64}()

        # Solve A: slab DL on bare steel — use section props from catalogue
        for i in 1:lastindex(beam_elements)
            j_sel = section_indices[i]
            beam_elements[i].section = Section(catalogue_df.A[j_sel], E_s, steel_ksi.G,
                catalogue_df.Ix[j_sel], catalogue_df.Iy[j_sel], catalogue_df.J[j_sel])
        end
        update_load_values_staged!(params.model, params, load_case=:slab_dead)
        params.load_dictionary = get_load_dictionary_by_id(params.model)
        Asap.solve!(params.model, reprocess=true)

        δ_dead_bare = Dict{Int, Float64}()
        bare_Ix_used = Dict{Int, Float64}()
        for i in 1:lastindex(beam_elements)
            beam_id = get_element_id(beam_elements[i])
            disp = ElementDisplacements(beam_elements[i], params.load_dictionary[beam_id], resolution=200)
            δ_dead_bare[i] = maximum(abs.(disp.ulocal[2, :]))
            bare_Ix_used[i] = beam_elements[i].section.Ix
        end

        # Solve B: SDL+LL on composite — use catalogue A/Iy/J with composite Ix
        for i in 1:lastindex(beam_elements)
            j_sel = section_indices[i]
            beam_elements[i].section = Section(catalogue_df.A[j_sel], E_s, steel_ksi.G,
                effective_Ix[i][j_sel], catalogue_df.Iy[j_sel], catalogue_df.J[j_sel])
        end
        update_load_values_staged!(params.model, params, load_case=:sdl_live)
        params.load_dictionary = get_load_dictionary_by_id(params.model)
        Asap.solve!(params.model, reprocess=true)

        δ_sdl_live_comp = Dict{Int, Float64}()
        comp_Ix_used = Dict{Int, Float64}()
        for i in 1:lastindex(beam_elements)
            beam_id = get_element_id(beam_elements[i])
            disp = ElementDisplacements(beam_elements[i], params.load_dictionary[beam_id], resolution=200)
            δ_sdl_live_comp[i] = maximum(abs.(disp.ulocal[2, :]))
            comp_Ix_used[i] = beam_elements[i].section.Ix
        end

        # Derive Ix requirements per beam using the same *proportional L/240
        # allocation* as the continuous NLP (see lines ~850 below). Splitting
        # lim_240 in proportion to each side's current deflection share makes
        # both bare and composite constraints scale by the same factor k =
        # δ_total / lim_240, so the MIP cannot trade bare against composite.
        #
        #     δ_bare_new ≤ lim_240 · f_bare,   δ_sdl_live_new ≤ lim_240 · f_comp
        #     f_bare = δ_bare_total_curr / δ_total_curr
        #     f_comp = δ_sdl_live_curr   / δ_total_curr   (= 1 − f_bare)
        #
        # L/360 (live only) is an independent composite constraint.
        loadid_index_stg = _build_loadid_index(params)
        for i in 1:lastindex(beam_elements)
            L_beam  = beam_elements[i].length
            lim_240 = L_beam / 240.0
            lim_360 = L_beam / 360.0

            w_sdl_sum = 0.0; w_live_sum = 0.0
            beam_id = get_element_id(beam_elements[i])
            for ld in params.load_dictionary[beam_id]
                if hasproperty(ld, :loadID)
                    row = get(loadid_index_stg, ld.loadID, nothing)
                    if !isnothing(row)
                        w_sdl_sum  += params.load_df[row, :unfactored_w_sdl]
                        w_live_sum += params.load_df[row, :unfactored_w_live]
                    end
                end
            end
            w_sum = w_sdl_sum + w_live_sum
            f_live_i = w_sum > 0 ? w_live_sum / w_sum : 0.5

            # L/360 live — independent, composite only
            δ_live_i = δ_sdl_live_comp[i] * f_live_i
            Ix_live_req[i] = δ_live_i > 0 ? comp_Ix_used[i] * δ_live_i / lim_360 : 0.0

            # Self-weight deflection at the current seed section
            j_sel = section_indices[i]
            δ_sw_i         = 5 * ρ_steel * catalogue_df.A[j_sel] * L_beam^4 /
                             (384 * E_s * max(bare_Ix_used[i], 1e-12))
            δ_bare_total_i = δ_dead_bare[i] + δ_sw_i
            δ_total_i      = δ_bare_total_i + δ_sdl_live_comp[i]

            # Proportional L/240 allocation (always active as a floor).
            if δ_total_i > 0 && δ_sdl_live_comp[i] > 0
                f_bare   = δ_bare_total_i / δ_total_i
                # Floor to keep either side from collapsing when one dominates
                # at a degenerate seed; 10 % gives the solver room to find an
                # interior allocation.
                f_bare   = clamp(f_bare, 0.1, 0.9)
                f_comp   = 1.0 - f_bare
                lim_bare = lim_240 * f_bare
                lim_comp = lim_240 * f_comp

                # Bare: (Ix_bare − sw_c·A) ≥ δ_dead_bare · bare_Ix_curr / lim_bare
                # (sw isolated so A can vary in the MIP section choice)
                Ix_coupled_bare[i] = δ_dead_bare[i] * bare_Ix_used[i] / lim_bare
                sw_coeff_240[i]    = 5 * ρ_steel * L_beam^4 /
                                     (384 * E_s * lim_bare)
                Ix_coupled_comp[i] = δ_sdl_live_comp[i] * comp_Ix_used[i] / lim_comp
            else
                # No composite-side load (e.g. no slab reaches this beam):
                # bare-only L/240 constraint on the dead-load stage.
                Ix_coupled_comp[i] = 0.0
                if δ_dead_bare[i] > 0
                    Ix_coupled_bare[i] = δ_dead_bare[i] * bare_Ix_used[i] / lim_240
                    sw_coeff_240[i]    = 5 * ρ_steel * L_beam^4 /
                                         (384 * E_s * lim_240)
                else
                    Ix_coupled_bare[i] = 0.0
                    sw_coeff_240[i]    = 0.0
                end
            end
        end

        # Restore factored loads
        update_load_values!(params.model, params, factored=true)
        params.load_dictionary = get_load_dictionary_by_id(params.model)
        Asap.solve!(params.model, reprocess=true)

        return Ix_live_req, Ix_coupled_comp, Ix_coupled_bare, sw_coeff_240
    end

    # ── legacy (non-staged) deflection requirements ───────────────────────
    Ix_base_req = Dict{Int, Float64}()
    sw_coeff_req = Dict{Int, Float64}()

    Ix_live_req       = Dict{Int, Float64}()
    Ix_coupled_comp   = Dict{Int, Float64}()
    Ix_coupled_bare   = Dict{Int, Float64}()
    sw_coeff_240      = Dict{Int, Float64}()

    if params.deflection_limit
        if use_staged_mip
            # Initial pass: use lightest feasible section for each beam
            init_section_indices = Dict(i => filtered_indices[i][1] for i in 1:lastindex(beam_elements))
            Ix_live_req, Ix_coupled_comp, Ix_coupled_bare, sw_coeff_240 =
                _compute_staged_Ix_requirements(init_section_indices)

            # Prefilter: L/360 (live) uses composite Ix; L/240 uses both sides
            # with a proportional split (composite via effective_Ix, bare via
            # catalogue Ix minus sw_c·A).
            tol = 1e-9
            for i in 1:lastindex(beam_elements)
                min_comp_ix = max(Ix_live_req[i], get(Ix_coupled_comp, i, 0.0))
                min_bare_ix = get(Ix_coupled_bare, i, 0.0)
                sw_c        = get(sw_coeff_240, i, 0.0)
                feasible = [j for j in filtered_indices[i] if
                    effective_Ix[i][j] + tol >= min_comp_ix &&
                    (catalogue_df.Ix[j] - sw_c * catalogue_df.A[j]) + tol >= min_bare_ix]
                if isempty(feasible)
                    throw(NoValidSectionsError("No sections satisfy staged deflection constraints for beam $i."))
                end
                filtered_indices[i] = feasible
            end
        else
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
                    heaviest = filtered_indices[i][argmax([effective_Ix[i][j] for j in filtered_indices[i]])]
                    @warn "Beam $i: no section satisfies staged Ix floor " *
                          "(need comp=$(round(floor_comp,digits=1)), bare=$(round(floor_bare,digits=1))); " *
                          "keeping heaviest feasible (comp=$(round(effective_Ix[i][heaviest],digits=1)))"
                    filtered_indices[i] = [heaviest]
                else
                    filtered_indices[i] = feasible
                end
            end
            max_comp_avail = maximum(effective_Ix[i][j] for j in filtered_indices[i])
            if use_staged_mip
                Ix_live_req[i]     = min(max(get(Ix_live_req, i, 0.0), floor_comp), max_comp_avail)
                Ix_coupled_comp[i] = min(max(get(Ix_coupled_comp, i, 0.0), floor_comp), max_comp_avail)
            elseif params.deflection_limit
                Ix_base_req[i] = min(max(get(Ix_base_req, i, 0.0), floor_comp), max_comp_avail)
            end
        end
    end

    # ── iterative MIP solve with staged deflection tightening ─────────────
    max_defl_iters = use_staged_mip ? 5 : 1
    chosen_section_idx = Dict{Int, Int}()

    for defl_iter in 1:max_defl_iters

        if defl_iter > 1
            println("  ── staged deflection re-solve (iteration $defl_iter) ──")
        end

        # Unify filtered indices within collinear groups
        groups = params.collinear_groups
        local filtered_iter = deepcopy(filtered_indices)
        if !isempty(groups)
            for gid in unique(groups)
                members = findall(==(gid), groups)
                if length(members) > 1
                    common = intersect([filtered_iter[m] for m in members]...)
                    if isempty(common)
                        @warn "Collinear group $gid has no common feasible sections — relaxing to union"
                        common = sort(unique(vcat([filtered_iter[m] for m in members]...)))
                    end
                    for m in members
                        filtered_iter[m] = common
                    end
                end
            end
        end

        jump_model = JuMP.Model(Gurobi.Optimizer)
        JuMP.set_optimizer_attribute(jump_model, "OutputFlag", defl_iter == 1 ? 1 : 0)
        JuMP.set_optimizer_attribute(jump_model, "MIPGap", 1e-4)
        JuMP.set_optimizer_attribute(jump_model, "NodeLimit", 10000)

        JuMP.@variable(jump_model, x[i=1:lastindex(beam_elements), j=filtered_iter[i]], Bin)
        JuMP.@constraint(jump_model, [i in 1:lastindex(beam_elements)],
            sum(x[i,j] for j in filtered_iter[i]) == 1)

        for i in 1:lastindex(beam_elements)
            JuMP.@constraint(jump_model, ϕ_b * sum(x[i,j] * catalogue_df.Mn[j] for j in filtered_iter[i]) >= beams_df.M_max[i])
            JuMP.@constraint(jump_model, ϕ_v * sum(x[i,j] * catalogue_df.Vn[j] for j in filtered_iter[i]) >= beams_df.V_max[i])
            JuMP.@constraint(jump_model, sum(x[i,j] * catalogue_df.h[j] for j in filtered_iter[i]) <= params.max_depth)
        end

        if params.deflection_limit
            if use_staged_mip
                for i in 1:lastindex(beam_elements)
                    JuMP.@constraint(jump_model,
                        sum(x[i,j] * effective_Ix[i][j] for j in filtered_iter[i]) >= Ix_live_req[i])
                    JuMP.@constraint(jump_model,
                        sum(x[i,j] * effective_Ix[i][j] for j in filtered_iter[i]) >= Ix_coupled_comp[i])
                    # Bare-side L/240 (proportional split): the dead-load stage
                    # must fit in its share of the span limit regardless of
                    # composite stiffness. Without this, the MIP can pick a
                    # section that passes the composite constraint but still
                    # violates L/240 on bare steel alone.
                    if get(Ix_coupled_bare, i, 0.0) > 0
                        sw_c = sw_coeff_240[i]
                        JuMP.@constraint(jump_model,
                            sum(x[i,j] * (catalogue_df.Ix[j] - sw_c * catalogue_df.A[j])
                                for j in filtered_iter[i]) >= Ix_coupled_bare[i])
                    end
                end
            else
                for i in 1:lastindex(beam_elements)
                    JuMP.@constraint(jump_model,
                        sum(x[i,j] * (effective_Ix[i][j] - sw_coeff_req[i] * catalogue_df.A[j])
                            for j in filtered_iter[i]) >= Ix_base_req[i])
                end
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

        if JuMP.termination_status(jump_model) != JuMP.MOI.OPTIMAL
            @warn "Optimization did not reach optimal solution. Status: $(JuMP.termination_status(jump_model))"
            return params
        end

        # Extract chosen section index per beam
        for i in 1:lastindex(beam_elements)
            chosen_section_idx[i] = filtered_iter[i][argmax([JuMP.value(x[i,j]) for j in filtered_iter[i]])]
        end

        # If not using staged iteration, one pass is enough
        if !use_staged_mip
            break
        end

        # Re-solve staged FE with chosen sections and check violations
        Ix_live_new, Ix_coupled_new, Ix_bare_new, sw_coeff_new =
            _compute_staged_Ix_requirements(chosen_section_idx)

        n_violations = 0
        tol = 1e-9
        for i in 1:lastindex(beam_elements)
            j_sel = chosen_section_idx[i]
            min_comp_needed = max(Ix_live_new[i], Ix_coupled_new[i])
            bare_margin = catalogue_df.Ix[j_sel] - sw_coeff_new[i] * catalogue_df.A[j_sel]
            if effective_Ix[i][j_sel] + tol < min_comp_needed ||
               bare_margin + tol < Ix_bare_new[i]
                n_violations += 1
            end
        end

        if n_violations == 0
            println("  Staged deflection converged after $defl_iter iteration(s) — no violations.")
            break
        end

        println("  $n_violations beam(s) violate staged deflection limits — updating requirements...")

        # Use latest requirements directly (not monotone max) so the MIP can
        # relax if the new FE solve shows the prior floor was too aggressive.
        for i in 1:lastindex(beam_elements)
            Ix_live_req[i]     = Ix_live_new[i]
            Ix_coupled_comp[i] = Ix_coupled_new[i]
            Ix_coupled_bare[i] = Ix_bare_new[i]
            sw_coeff_240[i]    = sw_coeff_new[i]
        end

        # Re-prefilter with updated requirements (both composite and bare).
        for i in 1:lastindex(beam_elements)
            min_comp_ix = max(Ix_live_req[i], Ix_coupled_comp[i])
            min_bare_ix = get(Ix_coupled_bare, i, 0.0)
            sw_c        = get(sw_coeff_240, i, 0.0)
            feasible = [j for j in filtered_indices[i] if
                effective_Ix[i][j] + tol >= min_comp_ix &&
                (catalogue_df.Ix[j] - sw_c * catalogue_df.A[j]) + tol >= min_bare_ix]
            if isempty(feasible)
                heaviest = filtered_indices[i][argmax([effective_Ix[i][j] for j in filtered_indices[i]])]
                @warn "Beam $i: no section satisfies staged constraints; keeping heaviest feasible"
                filtered_indices[i] = [heaviest]
            else
                filtered_indices[i] = feasible
            end
        end

        if defl_iter == max_defl_iters
            @warn "Staged deflection did not fully converge after $max_defl_iters iterations; $n_violations violation(s) remain."
        end
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

function process_continuous_beams_parallel(params::SlabSizingParams;
                                           initial_vars::Vector=[],
                                           max_demand_iters::Int=5,
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
        local Ix_req_live_vec, Ix_req_coupled_vec, Ix_req_bare_vec, sw_coeff_vec
        local use_staged_cont
        if params.deflection_limit
            has_staged_cont = :unfactored_w_live in propertynames(params.load_df)
            use_staged_cont = has_staged_cont && params.composite_action

            if use_staged_cont
                ρ_steel_cont = steel_ksi.ρ
                E_s = steel_ksi.E

                # Compute fresh section properties from current_vars
                fresh_A  = Vector{Float64}(undef, n_beams)
                fresh_Ix = Vector{Float64}(undef, n_beams)
                fresh_Iy = Vector{Float64}(undef, n_beams)
                fresh_J  = Vector{Float64}(undef, n_beams)
                for i in 1:n_beams
                    hi, wi, twi, tfi = current_vars[i]
                    fresh_A[i]  = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
                    fresh_Ix[i] = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
                    fresh_Iy[i] = Iy_I_asymm(hi, wi, wi, twi, tfi, tfi)
                    fresh_J[i]  = J_I_symm(hi, wi, twi, tfi)
                end

                # Solve A: slab DL on bare steel — set sections from current_vars
                for i in 1:n_beams
                    beam_elements[i].section = Section(fresh_A[i], E_s, steel_ksi.G,
                        fresh_Ix[i], fresh_Iy[i], fresh_J[i])
                end
                update_load_values_staged!(params.model, params, load_case=:slab_dead)
                params.load_dictionary = get_load_dictionary_by_id(params.model)
                Asap.solve!(params.model, reprocess=true)

                δ_dead_bare = Vector{Float64}(undef, n_beams)
                bare_Ix_used = Vector{Float64}(undef, n_beams)
                for i in 1:n_beams
                    beam_id = get_element_id(beam_elements[i])
                    disp = ElementDisplacements(beam_elements[i], params.load_dictionary[beam_id], resolution=200)
                    δ_dead_bare[i] = maximum(abs.(disp.ulocal[2, :]))
                    bare_Ix_used[i] = beam_elements[i].section.Ix
                end

                # Solve B: SDL + LL on composite — use fresh A/Iy/J with composite Ix
                for i in 1:n_beams
                    hi, wi, twi, tfi = current_vars[i]
                    Ix_comp_i = fresh_Ix[i]
                    if params.composite_action && params.slab_depth_in > 0
                        positions_i, widths_i, is_perim, L_beam = beam_load_data_cont[i]
                        if !isempty(positions_i)
                            Ix_comp_i = get_I_composite_effective(hi, wi, twi, tfi,
                                params.slab_depth_in, E_s, params.E_c, L_beam,
                                positions_i, widths_i; is_perimeter=is_perim)
                        end
                    end
                    beam_elements[i].section = Section(fresh_A[i], E_s, steel_ksi.G,
                        Ix_comp_i, fresh_Iy[i], fresh_J[i])
                end
                update_load_values_staged!(params.model, params, load_case=:sdl_live)
                params.load_dictionary = get_load_dictionary_by_id(params.model)
                Asap.solve!(params.model, reprocess=true)

                δ_sdl_live = Vector{Float64}(undef, n_beams)
                comp_Ix_used = Vector{Float64}(undef, n_beams)
                for i in 1:n_beams
                    beam_id = get_element_id(beam_elements[i])
                    disp = ElementDisplacements(beam_elements[i], params.load_dictionary[beam_id], resolution=200)
                    δ_sdl_live[i] = maximum(abs.(disp.ulocal[2, :]))
                    comp_Ix_used[i] = beam_elements[i].section.Ix
                end

                # Restore bare-steel sections and factored loads
                for i in 1:n_beams
                    beam_elements[i].section = Section(fresh_A[i], E_s, steel_ksi.G,
                        fresh_Ix[i], fresh_Iy[i], fresh_J[i])
                end
                update_load_values!(params.model, params, factored=true)
                params.load_dictionary = get_load_dictionary_by_id(params.model)
                Asap.solve!(params.model, reprocess=true)

                # Derive Ix requirements per beam using a *proportional L/240
                # allocation*. The one-sided residual-budget approach the MIP
                # uses is not self-consistent when the NLP is free to shrink
                # bare Ix: dead-load deflection grows and eats the composite
                # side's budget, but the NLP sees a stale residual and declares
                # feasibility. Instead, split lim_240 between bare and comp
                # *proportionally* to each side's current-section contribution:
                #
                #     δ_bare_total_new ≤ lim_240 · f_bare,   δ_sdl_live_new ≤ lim_240 · f_comp
                #     f_bare = δ_bare_total_curr / δ_total_curr
                #     f_comp = δ_sdl_live_curr  / δ_total_curr  (= 1 − f_bare)
                #
                # Each Ix must scale by the same factor k = δ_total / lim_240, so
                # the NLP cannot trade bare against composite. When k ≤ 1 the
                # section already satisfies L/240 and no tightening is needed.
                Ix_req_live_vec    = Vector{Float64}(undef, n_beams)
                Ix_req_coupled_vec = Vector{Float64}(undef, n_beams)
                Ix_req_bare_vec    = Vector{Float64}(undef, n_beams)
                sw_coeff_vec       = Vector{Float64}(undef, n_beams)
                loadid_index_defl  = _build_loadid_index(params)

                for i in 1:n_beams
                    L_beam  = beam_elements[i].length
                    lim_240 = L_beam / 240.0
                    lim_360 = L_beam / 360.0

                    # Self-weight deflection at the current bare section
                    δ_sw_i         = 5 * ρ_steel_cont * fresh_A[i] * L_beam^4 /
                                     (384 * E_s * max(bare_Ix_used[i], 1e-12))
                    δ_bare_total_i = δ_dead_bare[i] + δ_sw_i
                    δ_total_i      = δ_bare_total_i + δ_sdl_live[i]

                    # Live-load fraction of the unfactored post-composite load
                    w_sdl_sum = 0.0; w_live_sum = 0.0
                    beam_id = get_element_id(beam_elements[i])
                    for ld in params.load_dictionary[beam_id]
                        if hasproperty(ld, :loadID)
                            row = get(loadid_index_defl, ld.loadID, nothing)
                            if !isnothing(row)
                                w_sdl_sum  += params.load_df[row, :unfactored_w_sdl]
                                w_live_sum += params.load_df[row, :unfactored_w_live]
                            end
                        end
                    end
                    w_sum    = w_sdl_sum + w_live_sum
                    f_live_i = w_sum > 0 ? w_live_sum / w_sum : 0.5

                    # L/360 live-only: composite Ix must cover the live-only part
                    δ_live_i = δ_sdl_live[i] * f_live_i
                    Ix_req_live_vec[i] = δ_live_i > 0 ?
                        comp_Ix_used[i] * δ_live_i / lim_360 : 0.0

                    # Proportional L/240 allocation is *always* active: it
                    # acts as a floor that scales with demand. When the
                    # current section is overkill (δ_total_curr < lim_240),
                    # the target Ix is below the current Ix, so the NLP can
                    # shrink down to the feasibility edge; when the section
                    # is tight, the target pushes both bare and composite up.
                    if δ_total_i > 0 && δ_sdl_live[i] > 0
                        f_bare   = δ_bare_total_i / δ_total_i
                        # Floor to keep either side from collapsing when one
                        # dominates at a degenerate seed; 10% gives the solver
                        # room to find an interior allocation.
                        f_bare   = clamp(f_bare, 0.1, 0.9)
                        f_comp   = 1.0 - f_bare
                        lim_bare = lim_240 * f_bare
                        lim_comp = lim_240 * f_comp

                        # Bare: Ix_bare − sw_c·A ≥ δ_dead_bare·bare_Ix_curr / lim_bare
                        # (sw isolated so A can vary at the new section)
                        Ix_req_bare_vec[i]    = δ_dead_bare[i] * bare_Ix_used[i] / lim_bare
                        sw_coeff_vec[i]       = 5 * ρ_steel_cont * L_beam^4 /
                                                (384 * E_s * lim_bare)
                        Ix_req_coupled_vec[i] = δ_sdl_live[i] * comp_Ix_used[i] / lim_comp
                    else
                        # No composite-side load (e.g. no slab): fall back to a
                        # pure bare-steel L/240 constraint on dead-only.
                        Ix_req_coupled_vec[i] = 0.0
                        if δ_dead_bare[i] > 0
                            Ix_req_bare_vec[i] = δ_dead_bare[i] * bare_Ix_used[i] / lim_240
                            sw_coeff_vec[i]    = 5 * ρ_steel_cont * L_beam^4 /
                                                 (384 * E_s * lim_240)
                        else
                            Ix_req_bare_vec[i] = 0.0
                            sw_coeff_vec[i]    = 0.0
                        end
                    end
                end
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
                Ix_req_live_vec    = zeros(n_beams)
                Ix_req_coupled_vec = zeros(n_beams)
                Ix_req_bare_vec    = Vector{Float64}(undef, n_beams)
                sw_coeff_vec       = Vector{Float64}(undef, n_beams)
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

        # Fix C: prevent optimizer from runaway growth past the warm-start.
        # The MIP warm-start is a feasible discrete solution; the NLP relaxes
        # to continuous sections and should generally improve or match it, but
        # demand-loop drift and numerical slack can briefly push the iterate
        # above. We cap at 10% growth to absorb this without hiding real
        # infeasibility — any larger excess signals the NLP is finding a
        # heavier "improvement" than the MIP, which is worth reviewing.
        if has_warm_start
            ws_obj = sum(warm_start_areas[i] * beams_df.x_max[i] for i in 1:n_beams)
            JuMP.@constraint(jump_model, sum(
                op_A_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) * beams_df.x_max[i]
                for i in 1:n_beams) <= ws_obj * 1.10)
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
        # Staged: three constraints per beam
        #   (a) L/360 live  — composite Ix ≥ Ix_req_live
        #   (b) L/240 total, composite share  — composite Ix ≥ Ix_req_coupled
        #   (c) L/240 total, bare share       — Ix_bare − sw_c·A ≥ Ix_req_bare
        #   (b) and (c) use a *proportional* allocation of lim_240; each side
        #   scales by the same factor k = δ_total/lim_240, so the NLP cannot
        #   undersize bare by overshooting composite (the failure mode that
        #   was producing L/240-violating "feasible" solutions).
        # Legacy: single Ix_req on effective Ix (composite via frozen ratio
        #         or bare), net of a self-weight deduction.
        if params.deflection_limit
            if use_staged_cont
                for i in 1:n_beams
                    sw_c    = sw_coeff_vec[i]
                    ratio_i = (needs_composite_ratio && composite_ratio_vec[i] > 1.0) ?
                              composite_ratio_vec[i] : 1.0
                    # Skip trivial constraints (RHS = 0) to reduce solver stress
                    # — especially helpful for MMA on small geometries where
                    # deflection demands are already zero (e.g., no slab load).
                    if Ix_req_live_vec[i] > 0
                        JuMP.@constraint(jump_model,
                            Ix_req_live_vec[i] -
                            ratio_i * op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
                    end
                    if Ix_req_coupled_vec[i] > 0
                        JuMP.@constraint(jump_model,
                            Ix_req_coupled_vec[i] -
                            ratio_i * op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i]) <= 0)
                    end
                    # Bare-side L/240 is only active when the current section
                    # exceeds the total limit (sw_c > 0 ⇔ Ix_req_bare > 0).
                    if sw_c > 0 && Ix_req_bare_vec[i] > 0
                        JuMP.@constraint(jump_model,
                            Ix_req_bare_vec[i] -
                            (op_Ix_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i])
                             - sw_c * op_A_I_asymm(h[i], w[i], w[i], tw[i], tf[i], tf[i])) <= 0)
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

        # Guard against silent failures: abort the demand loop if the solver
        # did not return a usable primal. LOCALLY_SOLVED and OPTIMAL are good;
        # anything else (INFEASIBLE, ITERATION_LIMIT, ALMOST_*…) means the
        # reported values may violate constraints. We fall back to the best
        # solution seen in prior iterations, or report no minimizers.
        term = JuMP.termination_status(jump_model)
        prim = JuMP.primal_status(jump_model)
        good = (term == JuMP.MOI.OPTIMAL || term == JuMP.MOI.LOCALLY_SOLVED) &&
               prim == JuMP.MOI.FEASIBLE_POINT
        if !good
            @warn "Continuous NLP returned non-optimal status" demand_iter term prim
            # Fall back to best_vars (initialized to the MIP warm-start
            # before the loop). Post-processing will flag any catastrophic
            # deflection outcomes via SERVICEABILITY_FAIL / result_ok=false,
            # so we don't need to mark the configuration infeasible here —
            # doing so would hide cases where MIP found a viable discrete
            # solution that the NLP simply couldn't refine further.
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

        # Fix B: keep the best mass-feasible solution seen so far. In staged
        # mode we only accept an iterate as "best" if it also satisfies the
        # proportional L/240 split implied by the constraint targets — i.e.
        # δ_total ≤ lim_240 × 1.02 at the reported section. This guards the
        # early-termination path on mass convergence from locking in an
        # iterate that still violates deflection.
        defl_ok = true
        if use_staged_cont
            for i in 1:n_beams
                hi, wi, twi, tfi = current_vars[i]
                A_new  = A_I_asymm(hi, wi, wi, twi, tfi, tfi)
                Ix_new = Ix_I_asymm(hi, wi, wi, twi, tfi, tfi)
                L_beam = beam_elements[i].length
                lim_240 = L_beam / 240.0
                δ_sw_new   = 5 * steel_ksi.ρ * A_new * L_beam^4 /
                             (384 * steel_ksi.E * max(Ix_new, 1e-12))
                # Linearized projection of staged deflection at the new section
                δ_dead_new = δ_dead_bare[i] * bare_Ix_used[i] / max(Ix_new, 1e-12)
                δ_live_new = δ_sdl_live[i]  * comp_Ix_used[i] / max(Ix_new, 1e-12) /
                             max(composite_ratio_vec[i], 1.0)
                δ_total_new = δ_dead_new + δ_sw_new + δ_live_new
                if δ_total_new > lim_240 * 1.02
                    defl_ok = false
                    break
                end
            end
        end

        if cur_obj < best_obj && defl_ok
            best_obj  = cur_obj
            best_vars = deepcopy(current_vars)
        end

        if demand_iter > 1 && rel_change < demand_tol && defl_ok
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

    for (i, load) in enumerate(beam_loads)
        # distributed load magnitudes in LCS
        px, py, pz = (R * load.value) .* [1, -1, -1]

        My .+= moment_function.(py, L, xinc, load.position)
        Vy .+= shear_function.(py, L, xinc, load.position)
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

    for (i, load) in enumerate(beam_loads)
        # distributed load magnitudes in LCS
        px, py, pz = (R * load.value) .* [1, -1, -1]

        My .+= moment_function.(py, L, xinc, load.position)
        Vy .+= shear_function.(py, L, xinc, load.position)
    end

    return abs.(My), abs.(Vy), collect(xinc)

end