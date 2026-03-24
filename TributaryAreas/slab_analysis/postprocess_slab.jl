"""
    postprocess_slab(self, params; check_collinear, resolution, concise)

Post-process beam sizing results: compute mass, embodied carbon, column sizing,
and staged deflection analysis.

# Staged deflection (unshored composite construction)

When `params.load_df` contains the decomposed load columns `unfactored_w_live`
and `unfactored_w_sdl` (added by `get_scaled_model`), **two** FE solves plus
one analytical term produce the full staged deflection breakdown:

| Component      | Method             | Stiffness       | Load case      |
|----------------|--------------------|-----------------|----------------|
| `δ_slab_dead`  | FE solve A         | Bare steel Ix   | Slab DL        |
| `δ_beam_dead`  | Analytic 5wL⁴/384EI| Bare steel Ix   | Self-weight    |
| `δ_sdl`        | FE solve B × ratio | Composite Ix    | SDL + LL split |
| `δ_live`       | FE solve B × ratio | Composite Ix    | SDL + LL split |

The combined SDL+LL solve is split by load-weight ratio (linear superposition).

Total deflection: `δ_total = δ_slab_dead + δ_beam_dead + δ_sdl + δ_live`

Deflection limits checked:
- **L/360** for live load only (`δ_live`)
- **L/240** for total load (`δ_total`)

Results are stored per beam in the returned `SlabOptimResults`.

# Arguments
- `self::SlabAnalysisParams`: Slab geometry and load data.
- `params::SlabSizingParams`: Beam sizing parameters (sections, load_df, model).
- `check_collinear`: Override collinear flag.
- `resolution`: Points per beam for deflection sampling (default 200).
- `concise`: Skip force/deflection reporting if true.

# Returns
- `SlabOptimResults` with mass, carbon, forces, column sizing, and staged deflections.
"""
function postprocess_slab(self::SlabAnalysisParams, params::SlabSizingParams; check_collinear::Union{Bool, Nothing}=nothing, resolution = 200, concise::Bool=false)
    
    beam_elements = params.model.elements[:beam]

    # Determine collinearity
    collinear = isnothing(check_collinear) ? params.collinear : check_collinear

    """# Check if all elements are the same type -- if not, set collinear to true
    element_types = unique(typeof.(beam_elements))
    if length(element_types) > 1
        collinear = true
    end"""

    # Collect collinear elements if necessary
    if collinear
        params.collinear_minimizers, params.collinear_ids, params.collinear_minimums = collect_collinear_elements(self, params)
        minimizers = params.collinear_minimizers
        ids = params.collinear_ids
        minimums = params.collinear_minimums
    else
        minimizers = params.minimizers
        ids = params.ids
        minimums = params.minimums
    end

    # Determine reinforcement ratio based on slab type
    ρ_reinforcement = self.slab_type in [:isotropic, :uniaxial] ? 0.01 : 0.02

    # Calculate volumes and masses
    if isempty(self.slab_depths) 
        println("Slab depths are empty")
    elseif isempty(self.areas)
        println("Slab areas are empty")
    end

    if length(self.slab_depths) == length(self.areas) # Determinate
        volume_slab_and_rebar = sum(self.slab_depths[i] * self.areas[i] for i in 1:lastindex(self.areas)) # m³
    else # Indeterminate
        volume_slab_and_rebar = self.area * self.slab_depths[1] # m³
    end
    volume_rebar = volume_slab_and_rebar * ρ_reinforcement # m³
    volume_slab = volume_slab_and_rebar - volume_rebar # m³

    ρ_conc = params.concrete_material.ρ_concrete     # kg/m³
    ecc_conc = params.concrete_material.ecc_concrete  # kgCO₂e/kg

    mass_slab = volume_slab * ρ_conc # kg
    norm_mass_slab = mass_slab / self.area # kg/m²
    embodied_carbon_slab = ecc_conc * norm_mass_slab

    mass_rebar = volume_rebar * ρ_REBAR * 1.05 # Add turnups and anchors
    norm_mass_rebar = mass_rebar / self.area # kg/m²
    embodied_carbon_rebar = ECC_REBAR * norm_mass_rebar

    n_beams = lastindex(beam_elements)
    volumes = Vector{Float64}(undef, n_beams)
    areas = Vector{Float64}(undef, n_beams)
    exposed_surface_areas = Vector{Float64}(undef, n_beams)
    xs = Vector{Vector{Float64}}(undef, n_beams)
    Ps = Vector{Vector{Float64}}(undef, n_beams)
    Pxs = Vector{Vector{Float64}}(undef, n_beams)
    Mys = Vector{Vector{Float64}}(undef, n_beams)
    Mns = Vector{Float64}(undef, n_beams)
    Vys = Vector{Vector{Float64}}(undef, n_beams)
    Vns = Vector{Float64}(undef, n_beams)

    beam_ids = [get_element_id(be) for be in beam_elements]

    # Ensure factored loads are active for internal force reporting and column sizing.
    # This is critical when postprocess_slab is called multiple times (e.g., noncollinear
    # then collinear) — the previous call may have left unfactored loads on the model.
    update_load_values!(params.model, params, factored=true)
    params.load_dictionary = get_load_dictionary_by_id(params.model)
    Asap.solve!(params.model, reprocess=true)

    for i in 1:n_beams
        beam_id = beam_ids[i]
        beam_loads = params.load_dictionary[beam_id]

        if !concise

            beam_forces = InternalForces(beam_elements[i], beam_loads, resolution=200)

            P, Px = Float64[], Float64[]
            x = beam_forces.x
            load_density = 100

            for load in beam_loads

                if is_pointload(load)
                    push!(P, load.value[3])
                    push!(Px, load.position * beam_elements[i].length)
                elseif is_lineload(load)
                    line_load_distributed = load.value[3] / load_density
                    positions = collect(0:1/(load_density-1):1)

                    for pos in positions
                        push!(P, line_load_distributed)
                        push!(Px, pos * beam_elements[i].length)
                    end
                end
            end

            My, Vy = beam_forces.My, beam_forces.Vy

            xs[i] = x
            Ps[i] = P
            Pxs[i] = Px
            Mys[i] = My
            Vys[i] = Vy

        end
            
        section = I_symm(minimizers[i]...)

        Mn = section.Mn
        Vn = section.Vn

        if beam_elements[i].nodeStart.nodeID == :wall && beam_elements[i].nodeEnd.nodeID == :wall
            section.A = 0.0
        end

        areas[i] = section.A
        
        volume = section.A * beam_elements[i].length

        Ix_eff = section.Ix
        if params.composite_action && params.slab_depth_in > 0
            L_beam = beam_elements[i].length
            element_loads_i = params.load_dictionary[beam_ids[i]]
            positions_i = Float64[]
            widths_i = Float64[]
            for ld in element_loads_i
                if hasproperty(ld, :loadID)
                    row = findfirst(==(getproperty(ld, :loadID)), params.load_df.loadID)
                    if !isnothing(row)
                        push!(positions_i, ld.position)
                        push!(widths_i, params.load_df[row, :trib_width])
                    end
                end
            end
            if !isempty(widths_i)
                is_perim = i in params.i_perimeter
                Ix_eff = get_I_composite_effective(section.h, section.w, section.tw, section.tf,
                                         params.slab_depth_in, steel_ksi.E, params.E_c,
                                         L_beam, positions_i, widths_i;
                                         is_perimeter=is_perim)
            end
        end

        beam_elements[i].section = Section(section.A, steel_ksi.E, steel_ksi.G, Ix_eff, section.Iy, section.J)
        exposed_surface_area = (2*section.h + 2*section.w - section.tw + 2*section.tf) * beam_elements[i].length 

        Mns[i] = Mn
        Vns[i] = Vn
        volumes[i] = volume
        exposed_surface_areas[i] = exposed_surface_area
    end

    # Calculate beam mass and embodied carbon
    mass_beams = sum(volumes) * convert_to_m[:in]^3 * ρ_STEEL
    println("mass of beams: ", round(mass_beams, digits=2), " kg")
    total_exposed_surface_area = sum(exposed_surface_areas) * convert_to_m[:in]^2
    println("total exposed surface area: ", round(total_exposed_surface_area, digits=2), "m²")
    fireproofing_volume = total_exposed_surface_area * 0.01
    println("fireproofing volume (10mm thick): ", round(fireproofing_volume, digits=2), " m³")
    fireproofing_mass = fireproofing_volume * 352 # kg/m³ https://www.isolatek.com/construction/commercial-products/medium-density/
    println("fireproofing mass (10mm thick): ", round(fireproofing_mass, digits=2), " kg")
    norm_mass_fireproofing = fireproofing_mass / self.area
    embodied_carbon_fireproofing = ECC_CONCRETE * norm_mass_fireproofing
    println("fireproofing normalized mass (10mm thick): ", round(norm_mass_fireproofing, digits=3), " kg/m²")
    println("fireproofing embodied carbon (10mm thick): ", round(embodied_carbon_fireproofing, digits=3), " kgCO₂eq/m²")
    norm_mass_beams = mass_beams / self.area
    println("slab surface area: ", round(self.area, digits=2), " m²")
    embodied_carbon_beams = ECC_STEEL * norm_mass_beams

    # --- Column sizing (uses factored loads, must run before unfactored re-solve) ---
    col_sections, col_Pu, col_ϕPn, col_util, mass_columns, norm_mass_columns, embodied_carbon_columns =
        size_columns(params.model, self.area)

    δ_locals = Vector{Vector{Float64}}(undef, n_beams)
    δ_globals = Vector{Vector{Float64}}(undef, n_beams)

    # =========================================================================
    # Staged deflection analysis (unshored composite construction)
    #
    # Physical model — consolidated into 2 FE solves + 1 analytical term:
    #   Solve A — bare-steel Ix  with slab DL loads
    #   Solve B — composite  Ix  with SDL + LL (combined), then split by ratio
    #   Analytic — beam self-weight δ = 5wL⁴/(384EI_bare), added to dead stage
    #
    # Deflection superposition:
    #   δ_total = δ_slab_dead + δ_beam_dead + δ_sdl + δ_live
    #
    # Limits (IBC / ASCE 7):
    #   L/360 for live load only  (δ_live)
    #   L/240 for total load      (δ_total)
    # =========================================================================

    δ_slab_dead_vec  = zeros(n_beams)
    δ_beam_dead_vec  = zeros(n_beams)
    δ_sdl_vec        = zeros(n_beams)
    δ_live_vec       = zeros(n_beams)
    δ_total_vec      = zeros(n_beams)
    Δ_lim_live_vec   = zeros(n_beams)
    Δ_lim_total_vec  = zeros(n_beams)
    δ_live_ok_vec    = fill(true, n_beams)
    δ_total_ok_vec   = fill(true, n_beams)

    composite_Ix = [beam_elements[i].section.Ix for i in 1:n_beams]

    # Precompute fresh section properties from minimizers for staged FE solves.
    # beam_elements may carry stale A/Iy/J from a previous call.
    fresh_secs = [I_symm(minimizers[i]...) for i in 1:n_beams]
    bare_Ix = [fresh_secs[i].Ix for i in 1:n_beams]
    sec_A   = [fresh_secs[i].A  for i in 1:n_beams]
    sec_Iy  = [fresh_secs[i].Iy for i in 1:n_beams]
    sec_J   = [fresh_secs[i].J  for i in 1:n_beams]

    has_staged_loads = :unfactored_w_live in propertynames(params.load_df)

    if has_staged_loads && !concise
        E_steel = steel_ksi.E
        ρ_steel_kip = steel_ksi.ρ

        # --- Solve A: Slab DL on bare steel ---
        for i in 1:n_beams
            beam_elements[i].section = Section(sec_A[i], E_steel, steel_ksi.G,
                bare_Ix[i], sec_Iy[i], sec_J[i])
        end
        update_load_values_staged!(params.model, params, load_case=:slab_dead)
        params.load_dictionary = get_load_dictionary_by_id(params.model)
        Asap.solve!(params.model, reprocess=true)

        for i in 1:n_beams
            disp = ElementDisplacements(beam_elements[i],
                params.load_dictionary[beam_ids[i]], resolution=resolution)
            δ_slab_dead_vec[i] = maximum(abs.(disp.ulocal[2, :]))
        end

        # --- Analytical beam self-weight deflection on bare steel ---
        # δ_beam = 5·w·L⁴ / (384·E·I_bare), simply-supported uniform load
        # w = A · ρ_steel  [kip/in],  L in inches
        for i in 1:n_beams
            w_sw = sec_A[i] * ρ_steel_kip
            L_in = beam_elements[i].length
            δ_beam_dead_vec[i] = 5 * w_sw * L_in^4 / (384 * E_steel * bare_Ix[i])
        end

        # --- Solve B: SDL + LL on composite section (single solve) ---
        for i in 1:n_beams
            beam_elements[i].section = Section(sec_A[i], E_steel, steel_ksi.G,
                composite_Ix[i], sec_Iy[i], sec_J[i])
        end
        update_load_values_staged!(params.model, params, load_case=:sdl_live)
        params.load_dictionary = get_load_dictionary_by_id(params.model)
        Asap.solve!(params.model, reprocess=true)

        δ_sdl_live_vec = zeros(n_beams)
        for i in 1:n_beams
            disp = ElementDisplacements(beam_elements[i],
                params.load_dictionary[beam_ids[i]], resolution=resolution)
            δ_sdl_live_vec[i] = maximum(abs.(disp.ulocal[2, :]))
        end

        # Split combined SDL+LL deflection by load ratio (linear elastic).
        # SDL and LL produce the same deflection shape (uniform on same beam),
        # so δ_sdl/δ_sdl_live = w_sdl/(w_sdl + w_live).
        loadid_index = _build_loadid_index(params)
        for i in 1:n_beams
            w_sdl_sum = 0.0
            w_live_sum = 0.0
            for ld in params.load_dictionary[beam_ids[i]]
                if hasproperty(ld, :loadID)
                    row = get(loadid_index, ld.loadID, nothing)
                    if !isnothing(row)
                        w_sdl_sum  += params.load_df[row, :unfactored_w_sdl]
                        w_live_sum += params.load_df[row, :unfactored_w_live]
                    end
                end
            end
            w_total = w_sdl_sum + w_live_sum
            frac_sdl = w_total > 0 ? w_sdl_sum / w_total : 0.5
            δ_sdl_vec[i]  = δ_sdl_live_vec[i] * frac_sdl
            δ_live_vec[i] = δ_sdl_live_vec[i] * (1.0 - frac_sdl)
        end

        # --- Superpose and check limits ---
        for i in 1:n_beams
            L = beam_elements[i].length
            δ_total_vec[i] = δ_slab_dead_vec[i] + δ_beam_dead_vec[i] +
                             δ_sdl_vec[i] + δ_live_vec[i]
            Δ_lim_live_vec[i]  = L / 360.0
            Δ_lim_total_vec[i] = L / 240.0
            δ_live_ok_vec[i]  = δ_live_vec[i] <= Δ_lim_live_vec[i]
            δ_total_ok_vec[i] = δ_total_vec[i] <= Δ_lim_total_vec[i]
        end

        n_live_fail  = count(.!δ_live_ok_vec)
        n_total_fail = count(.!δ_total_ok_vec)
        println("Staged deflection: $(n_live_fail) beams exceed L/360 (live), " *
                "$(n_total_fail) beams exceed L/240 (total)")
    end

    # --- Legacy total-unfactored deflection (Δ_local, Δ_global) ---
    # Use fresh A/Iy/J from minimizers with composite Ix, unfactored loads.
    for i in 1:n_beams
        beam_elements[i].section = Section(sec_A[i], steel_ksi.E, steel_ksi.G,
            composite_Ix[i], sec_Iy[i], sec_J[i])
    end
    update_load_values!(params.model, params, factored=false)
    params.load_dictionary = get_load_dictionary_by_id(params.model)
    Asap.solve!(params.model, reprocess=true)

    if !concise
        for i in 1:n_beams
            displacements = ElementDisplacements(beam_elements[i],
                params.load_dictionary[beam_ids[i]], resolution=resolution)
            δ_locals[i] = displacements.ulocal[2, :]
            δ_globals[i] = displacements.uglobal[2, :]
        end
    end

    # Restore factored loads so subsequent calls (e.g., collinear postprocessing,
    # column sizing in a second invocation) see the correct load state.
    update_load_values!(params.model, params, factored=true)
    params.load_dictionary = get_load_dictionary_by_id(params.model)
    Asap.solve!(params.model, reprocess=true)

    # --- Aggregate utilization metrics ---
    _max_util_M = 0.0
    _max_util_V = 0.0
    if !concise
        for i in 1:n_beams
            if Mns[i] > 0 && !isempty(Mys[i])
                _max_util_M = max(_max_util_M, maximum(abs.(Mys[i])) / Mns[i])
            end
            if Vns[i] > 0 && !isempty(Vys[i])
                _max_util_V = max(_max_util_V, maximum(abs.(Vys[i])) / Vns[i])
            end
        end
    end
    _max_col_util = isempty(col_util) ? 0.0 : maximum(col_util)
    _i_L360_fail  = findall(.!δ_live_ok_vec)
    _i_L240_fail  = findall(.!δ_total_ok_vec)
    _n_L360_fail  = length(_i_L360_fail)
    _n_L240_fail  = length(_i_L240_fail)

    _max_δ_total  = isempty(δ_total_vec) ? 0.0 : maximum(δ_total_vec)
    _max_bay_span = params.max_bay_span
    _global_δ_ok  = _max_bay_span <= 0 || _max_δ_total <= _max_bay_span / 180.0

    _nlp_solver_str = params.beam_sizer == :discrete ? "MIP" : String(params.nlp_solver)

    results = SlabOptimResults(
        slab_name              = self.slab_name,
        slab_type              = self.slab_type,
        vector_1d              = self.vector_1d,
        slab_sizer             = self.slab_sizer,
        beam_sizer             = params.beam_sizer,
        area                   = self.area,
        minimizers             = minimizers,
        minimums               = minimums,
        ids                    = ids,
        collinear              = collinear,
        max_depth              = params.max_depth,
        mass_beams             = mass_beams,
        norm_mass_beams        = norm_mass_beams,
        embodied_carbon_beams  = embodied_carbon_beams,
        mass_slab              = mass_slab,
        norm_mass_slab         = norm_mass_slab,
        embodied_carbon_slab   = embodied_carbon_slab,
        mass_rebar             = mass_rebar,
        norm_mass_rebar        = norm_mass_rebar,
        embodied_carbon_rebar  = embodied_carbon_rebar,
        areas                  = collect(Float64, areas),
        x                      = xs,
        P                      = Ps,
        Px                     = Pxs,
        My                     = Mys,
        Mn                     = collect(Float64, Mns),
        Vy                     = Vys,
        Vn                     = collect(Float64, Vns),
        Δ_local                = δ_locals,
        Δ_global               = δ_globals,
        sections               = collect(String, ids),
        col_sections           = col_sections,
        col_Pu                 = col_Pu,
        col_ϕPn                = col_ϕPn,
        col_util               = col_util,
        mass_columns           = mass_columns,
        norm_mass_columns      = norm_mass_columns,
        embodied_carbon_columns= embodied_carbon_columns,
        mass_fireproofing      = fireproofing_mass,
        norm_mass_fireproofing = norm_mass_fireproofing,
        embodied_carbon_fireproofing = embodied_carbon_fireproofing,
        δ_slab_dead            = δ_slab_dead_vec,
        δ_beam_dead            = δ_beam_dead_vec,
        δ_sdl                  = δ_sdl_vec,
        δ_live                 = δ_live_vec,
        δ_total                = δ_total_vec,
        Δ_limit_live           = Δ_lim_live_vec,
        Δ_limit_total          = Δ_lim_total_vec,
        δ_live_ok              = δ_live_ok_vec,
        δ_total_ok             = δ_total_ok_vec,
        max_δ_total            = _max_δ_total,
        max_bay_span           = _max_bay_span,
        global_δ_ok            = _global_δ_ok,
        max_util_M             = _max_util_M,
        max_util_V             = _max_util_V,
        max_col_util           = _max_col_util,
        n_L360_fail            = _n_L360_fail,
        i_L360_fail            = _i_L360_fail,
        n_L240_fail            = _n_L240_fail,
        i_L240_fail            = _i_L240_fail,
        nlp_solver             = _nlp_solver_str,
        deflection_limit       = params.deflection_limit,
        composite_action       = params.composite_action,
        staged_converged       = params.staged_converged,
        staged_n_violations    = params.staged_n_violations,
    )

    return results
end

"""
    size_columns(model, floor_area; storey_height_m, Kx, Ky, Lx_m, Ly_m)

Size each column individually for gravity-only axial compression (AISC 360 Ch. E).
Uses the factored axial force already in the solved model, searches the W-shape
catalog lightest-first, and returns vectors of results.

# Effective length parameters
- `Kx`, `Ky` : effective length factors (strong / weak axis, default 1.0)
- `Lx_m`, `Ly_m` : unbraced lengths per axis [m].  Both default to `storey_height_m`
  so a single-storey column uses the same length on both axes unless overridden.

Assumes the model is solved with **factored** loads at the time of call (the column
sizing runs before the unfactored-load serviceability re-solve).
"""
function size_columns(model::Asap.Model, floor_area::Float64;
                      storey_height_m::Float64=4.0,
                      Kx::Float64=1.0, Ky::Float64=1.0,
                      Lx_m::Float64=storey_height_m,
                      Ly_m::Float64=storey_height_m)
    col_elements = try model.elements[:column] catch; Element[] end
    n_cols = length(col_elements)

    if n_cols == 0
        return String[], Float64[], Float64[], Float64[], 0.0, 0.0, 0.0
    end

    Lx_in = Lx_m / convert_to_m[:in]
    Ly_in = Ly_m / convert_to_m[:in]
    storey_height_in = storey_height_m / convert_to_m[:in]

    catalog = allW_imperial()

    col_sections = Vector{String}(undef, n_cols)
    col_Pu      = Vector{Float64}(undef, n_cols)
    col_ϕPn_out = Vector{Float64}(undef, n_cols)
    col_util    = Vector{Float64}(undef, n_cols)
    col_vol_in3 = Vector{Float64}(undef, n_cols)

    for i in 1:n_cols
        element = col_elements[i]
        forces = InternalForces(element, Asap.AbstractLoad[], resolution=2)
        Pu_kip = maximum(abs.(forces.P))

        col_Pu[i] = Pu_kip

        found = false
        for w in catalog
            sec = I_symm(w.d, w.bf, w.tw, w.tf)
            ϕPn = get_ϕPn(sec, Lx_in, Ly_in; Kx=Kx, Ky=Ky)
            if ϕPn >= Pu_kip
                col_sections[i] = w.name
                col_ϕPn_out[i]  = ϕPn
                col_util[i]     = Pu_kip / ϕPn
                col_vol_in3[i]  = sec.A * storey_height_in
                found = true
                break
            end
        end

        if !found
            heaviest = catalog[end]
            sec = I_symm(heaviest.d, heaviest.bf, heaviest.tw, heaviest.tf)
            ϕPn = get_ϕPn(sec, Lx_in, Ly_in; Kx=Kx, Ky=Ky)
            col_sections[i] = heaviest.name * " (OVERSTRESSED)"
            col_ϕPn_out[i]  = ϕPn
            col_util[i]     = Pu_kip / ϕPn
            col_vol_in3[i]  = sec.A * storey_height_in
        end

        println("  Column $i: $(col_sections[i])  Pu=$(round(col_Pu[i], digits=1)) kip  " *
                "ϕPn=$(round(col_ϕPn_out[i], digits=1)) kip  util=$(round(col_util[i]*100, digits=1))%")
    end

    total_vol_m3 = sum(col_vol_in3) * convert_to_m[:in]^3
    mass_columns = total_vol_m3 * ρ_STEEL
    norm_mass_columns = mass_columns / floor_area
    embodied_carbon_columns = ECC_STEEL * norm_mass_columns

    println("Column steel: $(round(norm_mass_columns, digits=2)) kg/m²  " *
            "EC: $(round(embodied_carbon_columns, digits=2)) kgCO₂eq/m²")

    return col_sections, col_Pu, col_ϕPn_out, col_util,
           mass_columns, norm_mass_columns, embodied_carbon_columns
end
