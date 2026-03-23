"""
Compare beam deflections across three stiffness assumptions for all topology geometries.

Architecture:
  1. Build geometry + analyze slab ONCE per JSON
  2. Size beams independently for each case (different composite/drf settings)
  3. Compute deflections from a fresh model copy for each case

Cases:
  1. Bare steel              — composite_action=false, drf=1.0
  2. Bare steel (2.5× factor)— composite_action=false, drf=2.5
  3. Composite action        — composite_action=true,  drf=1.0
"""

include("_scripts.jl")
using PrettyTables

CairoMakie.activate!()

# ── colour palette ────────────────────────────────────────────────────────────
色 = Dict(
    :powderblue   => colorant"#aeddf5",
    :skyblue      => colorant"#70cbfd",
    :gold         => colorant"#df7e00",
    :magenta      => colorant"#dc267f",
    :orange       => colorant"#e04600",
    :ceruleanblue => colorant"#00AEEF",
    :charcoalgrey => colorant"#3e3e3e",
    :irispurple   => colorant"#4c2563",
    :darkpurple   => colorant"#130039",
    :lilac        => colorant"#A678B5",
)

cases = [
    (label="Bare steel",        color=色[:charcoalgrey], composite=false, drf=1.0),
    (label="Bare steel (2.5×)", color=色[:lilac],        composite=false, drf=2.5),
    (label="Composite",         color=色[:ceruleanblue], composite=true,  drf=1.0),
]

fontsize      = 11
smallfontsize = 8

# ── paths ─────────────────────────────────────────────────────────────────────
main_path = "Geometries/topology/"
save_path = "SlabDesignFactors/plot/figures/composite/"
mkpath(save_path)

json_files = filter(x -> endswith(x, ".json"), readdir(main_path))
isempty(json_files) && @warn "No JSON files found in $main_path"

# ── helpers ───────────────────────────────────────────────────────────────────

function order_group_along_chain(indices::Vector{Int}, elements)
    isempty(indices) && return Tuple{Int,Bool}[]
    ordered = Tuple{Int,Bool}[(indices[1], true)] # (element index, forward along element local x)
    remaining = setdiff(Set(indices), Set(first.(ordered)))
    first_id = elements[ordered[1][1]].nodeStart.nodeID
    last_id  = elements[ordered[1][1]].nodeEnd.nodeID
    while !isempty(remaining)
        found = false
        for k in remaining
            e = elements[k]
            if e.nodeStart.nodeID == last_id
                push!(ordered, (k, true)); last_id = e.nodeEnd.nodeID; found = true
            elseif e.nodeEnd.nodeID == first_id
                pushfirst!(ordered, (k, true)); first_id = e.nodeStart.nodeID; found = true
            elseif e.nodeEnd.nodeID == last_id
                push!(ordered, (k, false)); last_id = e.nodeStart.nodeID; found = true
            elseif e.nodeStart.nodeID == first_id
                pushfirst!(ordered, (k, false)); first_id = e.nodeEnd.nodeID; found = true
            end
            if found; pop!(remaining, k); break; end
        end
        found || break
    end
    return ordered
end

"""Chain beams in a group, normalize x to 0→1, subtract linear baseline."""
function group_deflection_curve(ordered_items::Vector{Tuple{Int,Bool}}, x_vecs, δ_vecs, beam_lengths)
    indices = first.(ordered_items)
    L_cum = vcat(0.0, cumsum([beam_lengths[i] for i in indices]))
    L_tot = L_cum[end]
    L_tot <= 0 && return (Float64[], Float64[])

    # Stitch each element in chain order and drop duplicate node points at interfaces.
    x_parts = Vector{Float64}[]
    δ_parts = Vector{Float64}[]
    for (j, (i, forward)) in enumerate(ordered_items)
        if forward
            x_local = x_vecs[i]
            δ_local = δ_vecs[i]
        else
            # Traverse this segment opposite to local element orientation.
            x_local = reverse(beam_lengths[i] .- x_vecs[i])
            δ_local = reverse(δ_vecs[i])
        end

        x_seg = L_cum[j] .+ x_local
        δ_seg = copy(δ_local)
        if j > 1
            x_seg = x_seg[2:end]
            δ_seg = δ_seg[2:end]
        end
        push!(x_parts, x_seg)
        push!(δ_parts, δ_seg)
    end

    x_global = vcat(x_parts...)
    δ_combined = vcat(δ_parts...)
    t = x_global ./ L_tot
    δ_start = δ_combined[1]
    δ_end   = δ_combined[end]
    δ_zeroed = δ_combined .- (δ_start .+ (δ_end - δ_start) .* t)
    return (t, δ_zeroed)
end

"""
Compute deflection vectors from a model given section minimizers.
Returns (x_vecs, δ_local, beam_lengths, mass_kg).

Switches to unfactored loads for serviceability and adds beam self-weight as
a LineLoad on each element. Uses composite Ix when `composite=true`.
"""
function compute_deflections(model, minimizers, params::SlabSizingParams;
                             composite::Bool=false, slab_depth_in::Float64=0.0,
                             E_c::Float64=57.0*sqrt(4000.0), resolution::Int=200)
    beam_elements = model.elements[:beam]
    load_dictionary = params.load_dictionary
    load_df = params.load_df
    beam_ids = [get_element_id(be) for be in beam_elements]
    n = length(beam_elements)

    ρ_steel = steel_ksi.ρ  # kip/in³

    volumes = Float64[]
    for i in 1:n
        section = I_symm(minimizers[i]...)

        if beam_elements[i].nodeStart.nodeID == :wall && beam_elements[i].nodeEnd.nodeID == :wall
            section.A = 0.0
        end

        Ix_eff = section.Ix
        if composite && slab_depth_in > 0
            L_beam = beam_elements[i].length
            element_loads_i = load_dictionary[beam_ids[i]]
            positions_i = Float64[]
            widths_i = Float64[]
            for ld in element_loads_i
                if hasproperty(ld, :loadID)
                    row = findfirst(==(getproperty(ld, :loadID)), load_df.loadID)
                    if !isnothing(row)
                        push!(positions_i, ld.position)
                        push!(widths_i, load_df[row, :trib_width])
                    end
                end
            end
            if !isempty(widths_i)
                is_perim = i in params.i_perimeter
                Ix_eff = get_I_composite_effective(section.h, section.w, section.tw, section.tf,
                                         slab_depth_in, steel_ksi.E, E_c,
                                         L_beam, positions_i, widths_i;
                                         is_perimeter=is_perim)
            end
        end

        beam_elements[i].section = Section(section.A, steel_ksi.E, steel_ksi.G,
                                           Ix_eff, section.Iy, section.J, ρ_steel)
        push!(volumes, section.A * beam_elements[i].length)
    end

    # Switch to unfactored loads for serviceability
    update_load_values!(model, params, factored=false)

    # Add beam self-weight as LineLoad on each element (unfactored dead load)
    sw_loads = Asap.AbstractLoad[]
    for be in beam_elements
        w_sw = be.section.A * ρ_steel  # kip/in (force per length)
        push!(sw_loads, LineLoad(be, [0.0, 0.0, -w_sw]))
    end
    append!(model.loads, sw_loads)

    Asap.solve!(model, reprocess=true)

    # Rebuild load dictionary to include self-weight loads
    load_dictionary_full = get_load_dictionary_by_id(model)

    x_vecs  = Vector{Vector{Float64}}(undef, n)
    δ_local = Vector{Vector{Float64}}(undef, n)
    for i in 1:n
        disp = ElementDisplacements(beam_elements[i], load_dictionary_full[beam_ids[i]],
                                    resolution=resolution)
        bf = InternalForces(beam_elements[i], load_dictionary_full[beam_ids[i]], resolution=resolution)
        x_vecs[i]  = bf.x
        # Use global vertical displacement so all elements in a collinear group
        # share one sign convention regardless of local element orientation.
        δ_local[i] = disp.uglobal[3, :]
    end

    beam_lengths = [beam_elements[i].length for i in 1:n]
    mass_kg = sum(volumes) * convert_to_m[:in]^3 * ρ_STEEL
    return x_vecs, δ_local, beam_lengths, mass_kg
end

# ── summary collection ────────────────────────────────────────────────────────
summary_df = DataFrame(
    geometry          = String[],
    n_beams           = Int[],
    total_beam_length = Float64[],
    max_δ_bare        = Float64[],
    max_δ_bare25      = Float64[],
    max_δ_composite   = Float64[],
    mass_bare         = Float64[],
    mass_bare25       = Float64[],
    mass_composite    = Float64[],
    # serviceability & strength (composite case)
    n_L360_fail       = Int[],
    n_L240_fail       = Int[],
    max_util_M        = Float64[],
    max_util_V        = Float64[],
    max_col_util      = Float64[],
    global_δ_ok       = Bool[],
)
all_max_δ = [Float64[] for _ in cases]
all_group_lengths = Float64[]
all_group_mass = [Float64[] for _ in cases]

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN LOOP
# ══════════════════════════════════════════════════════════════════════════════
for json_file in json_files

    name = replace(json_file, ".json" => "")
    path = joinpath(main_path, json_file)

    println("\n════════════════════════════════════════════")
    println("  $name")
    println("════════════════════════════════════════════")

    raw = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))
    geometry_dict = raw isa Dict ? raw : Dict(pairs(raw))

    # ── 1. Build geometry + analyze slab ONCE ─────────────────────────────
    geom, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)
    beam_elements_ref = geom.elements[:beam]
    n_beam = length(beam_elements_ref)
    beam_lengths = [be.length for be in beam_elements_ref]

    # Collinear groups (topology is identical across cases)
    collinear_groups = get_collinear_groups(beam_elements_ref)

    group_ids = unique(collinear_groups)
    group_ordered = [order_group_along_chain(findall(==(g), collinear_groups), beam_elements_ref) for g in group_ids]

    # ── 2. Size beams + postprocess for each case ─────────────────────────
    case_minimizers = Vector{Vector{Vector{Float64}}}(undef, length(cases))
    case_masses     = Vector{Float64}(undef, length(cases))
    case_results    = Vector{SlabOptimResults}(undef, length(cases))
    case_slab_depth = 0.0  # populated by composite case

    for (ci, c) in enumerate(cases)
        geom_c, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)

        slab_params = SlabAnalysisParams(
            geom_c,
            slab_name       = name,
            slab_type       = :isotropic,
            vector_1d       = [1.0, 0.0],
            slab_sizer      = :uniform,
            spacing         = 0.1,
            plot_analysis   = false,
            fix_param       = true,
            slab_units      = :m,
        )

        sizing_params = SlabSizingParams(
            live_load                   = psf_to_ksi(50),
            superimposed_dead_load      = psf_to_ksi(15),
            slab_dead_load              = 0.0,
            live_factor                 = 1.6,
            dead_factor                 = 1.2,
            beam_sizer                  = :discrete,
            max_depth                   = 40.0,
            beam_units                  = :in,
            serviceability_lim          = 360,
            collinear                   = true,
            minimum_continuous          = true,
            n_max_sections              = 0,
            composite_action            = c.composite,
            deflection_reduction_factor = c.drf,
        )

        slab_params = analyze_slab(slab_params)

        t0 = time()
        slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)
        dt = time() - t0

        case_minimizers[ci] = sizing_params.minimizers

        if isempty(sizing_params.minimizers)
            case_masses[ci] = 0.0
            case_results[ci] = SlabOptimResults()
        else
            scaled_beam_lengths = [be.length for be in sizing_params.model.elements[:beam]]
            case_masses[ci] = sum(
                I_symm(m...).A * scaled_beam_lengths[i] for (i, m) in enumerate(sizing_params.minimizers)
            ) * convert_to_m[:in]^3 * ρ_STEEL
            case_results[ci] = postprocess_slab(slab_params, sizing_params, check_collinear=true)
        end

        if c.composite
            case_slab_depth = sizing_params.slab_depth_in
        end

        println("  [$(c.label)] sized in $(round(dt, digits=2))s  " *
                "mass=$(round(case_masses[ci], digits=1))kg")
    end

    # ── 3. Compute deflections from fresh model copies ────────────────────
    # We need a model in imperial (inches) with the correct load dictionary.
    # Easiest: build a fresh model, run through the same scaling as optimal_beamsizer,
    # then call compute_deflections with each set of minimizers.

    case_x      = Vector{Vector{Vector{Float64}}}(undef, length(cases))
    case_δ      = Vector{Vector{Vector{Float64}}}(undef, length(cases))
    case_bL     = Vector{Vector{Float64}}(undef, length(cases))
    case_δ_mass = Vector{Float64}(undef, length(cases))

    for (ci, c) in enumerate(cases)
        geom_d, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)

        slab_params_d = SlabAnalysisParams(
            geom_d,
            slab_name       = name,
            slab_type       = :isotropic,
            vector_1d       = [1.0, 0.0],
            slab_sizer      = :uniform,
            spacing         = 0.1,
            plot_analysis   = false,
            fix_param       = true,
            slab_units      = :m,
        )

        sizing_params_d = SlabSizingParams(
            live_load                   = psf_to_ksi(50),
            superimposed_dead_load      = psf_to_ksi(15),
            slab_dead_load              = 0.0,
            live_factor                 = 1.6,
            dead_factor                 = 1.2,
            beam_sizer                  = :discrete,
            max_depth                   = 40.0,
            beam_units                  = :in,
            serviceability_lim          = 360,
            collinear                   = true,
            minimum_continuous          = true,
            n_max_sections              = 0,
            composite_action            = c.composite,
            deflection_reduction_factor = c.drf,
        )

        slab_params_d = analyze_slab(slab_params_d)

        # Scale model to inches (same transform as optimal_beamsizer)
        conversion_factor = convert_to_m[slab_params_d.slab_units] * 1/convert_to_m[sizing_params_d.beam_units]
        sizing_params_d.area = slab_params_d.area * conversion_factor^2
        sizing_params_d.model = get_scaled_model(slab_params_d, sizing_params_d, conversion_factor)
        sizing_params_d.load_dictionary = get_load_dictionary_by_id(sizing_params_d.model)

        if c.composite
            sizing_params_d.slab_depth_in = maximum(slab_params_d.slab_depths) * conversion_factor
            sizing_params_d.i_perimeter = Set(slab_params_d.i_perimeter)
        end

        x, δ, bL, mass = compute_deflections(
            sizing_params_d.model, case_minimizers[ci], sizing_params_d;
            composite   = c.composite,
            slab_depth_in = c.composite ? sizing_params_d.slab_depth_in : 0.0,
            E_c         = Float64(sizing_params_d.E_c),
        )

        case_x[ci]  = x
        case_δ[ci]  = δ
        case_bL[ci] = bL
        case_δ_mass[ci] = mass
    end

    # ── 4. Plot deflection profiles ───────────────────────────────────────
    fig = Figure(size=(800, 400), fontsize=fontsize)
    ax = Axis(fig[1, 1],
        xlabel = "Normalized length x/L",
        ylabel = "Deflection δ (in)",
        xlabelsize = smallfontsize, ylabelsize = smallfontsize,
        xticklabelsize = smallfontsize, yticklabelsize = smallfontsize,
    )

    δ_min_all =  Inf
    δ_max_all = -Inf
    for (ci, c) in enumerate(cases)
        for ordered in group_ordered
            xn, δ = group_deflection_curve(ordered, case_x[ci], case_δ[ci], case_bL[ci])
            isempty(xn) && continue
            δ_plot = ci == 2 ? δ ./ c.drf : δ
            lines!(ax, xn, δ_plot; color=c.color, linewidth=1.0,
                   linestyle=:solid)
            δ_min_all = min(δ_min_all, minimum(δ_plot))
            δ_max_all = max(δ_max_all, maximum(δ_plot))
        end
    end

    xlims!(ax, (0, 1))
    ylims!(ax, (δ_min_all * 1.1, δ_max_all * 1.1))
    for (ci, c) in enumerate(cases)
        lab = ci == 2 ? "$(c.label) / $(c.drf)" : c.label
        lines!(ax, [NaN], [NaN]; color=c.color, linewidth=2,
               linestyle=:solid, label=lab)
    end
    Legend(fig[1, 2], ax; labelsize=smallfontsize, framevisible=true,
           backgroundcolor=(:white, 0.9), framecolor=:grey80)
    Label(fig[0, :], "$name — deflection profiles by group",
          fontsize=fontsize, font=:bold, halign=:center)
    save(joinpath(save_path, "$(name)_deflection_profiles.pdf"), fig)

    # ── 5. Per-group max δ and per-group mass ───────────────────────────
    geom_max_δ = [0.0 for _ in cases]
    for ordered in group_ordered
        indices = first.(ordered)
        g_len = sum(case_bL[1][indices])
        push!(all_group_lengths, g_len)
        for ci in eachindex(cases)
            _, δ_zeroed = group_deflection_curve(ordered, case_x[ci], case_δ[ci], case_bL[ci])
            m = isempty(δ_zeroed) ? 0.0 : maximum(abs.(δ_zeroed))
            push!(all_max_δ[ci], m)
            geom_max_δ[ci] = max(geom_max_δ[ci], m)
            # Per-group mass (kg): sum of beam volumes × steel density
            g_mass = sum(I_symm(case_minimizers[ci][i]...).A * case_bL[ci][i] for i in indices) *
                     convert_to_m[:in]^3 * ρ_STEEL
            push!(all_group_mass[ci], g_mass)
        end
    end

    rc = case_results[3]  # composite case
    push!(summary_df, (
        name, n_beam, sum(beam_lengths),
        geom_max_δ[1], geom_max_δ[2], geom_max_δ[3],
        case_masses[1], case_masses[2], case_masses[3],
        rc.n_L360_fail, rc.n_L240_fail,
        rc.max_util_M, rc.max_util_V, rc.max_col_util,
        rc.global_δ_ok,
    ))

    for (ci, c) in enumerate(cases)
        println("  [$(c.label)] max δ: $(round(geom_max_δ[ci], digits=4)) in  " *
                "mass: $(round(case_masses[ci], digits=1)) kg")
    end
    if rc.n_L360_fail > 0
        println("  ⚠ L/360 fail ($(rc.n_L360_fail) beams): indices $(rc.i_L360_fail)")
    end
    if rc.n_L240_fail > 0
        println("  ⚠ L/240 fail ($(rc.n_L240_fail) beams): indices $(rc.i_L240_fail)")
    end
    if !rc.global_δ_ok
        println("  ⚠ Global δ sanity: max_δ_total=$(round(rc.max_δ_total, digits=4)) > span/180=$(round(rc.max_bay_span/180, digits=4))")
    end

    GC.gc()
end

# ══════════════════════════════════════════════════════════════════════════════
#  POST-LOOP: SUMMARY TABLE + PLOTS
# ══════════════════════════════════════════════════════════════════════════════

println("\n\n╔══════════════════════════════════════════════════════════════╗")
println("║           COMPOSITE ACTION COMPARISON SUMMARY              ║")
println("╚══════════════════════════════════════════════════════════════╝\n")

pretty_table(summary_df,
    column_labels = ["Geometry", "# Beams", "Total L",
                     "max δ bare", "max δ 2.5×", "max δ comp",
                     "Mass bare", "Mass 2.5×", "Mass comp",
                     "L/360 fail", "L/240 fail",
                     "max M util", "max V util", "max col util",
                     "global δ ok"],
    maximum_number_of_columns = -1,
    maximum_number_of_rows = -1,
    formatters = (v, i, j) -> begin
        if j in (4, 5, 6)
            return round(v, digits=4)
        elseif j in (7, 8, 9)
            return round(v, digits=1)
        elseif j in (12, 13, 14)
            return round(v, digits=3)
        else
            return v
        end
    end,
)

# ── summary bar chart ─────────────────────────────────────────────────────────
if nrow(summary_df) > 0
    ng = nrow(summary_df)
    fig_summary = Figure(size=(max(700, 90 * ng), 550), fontsize=fontsize)
    offsets = [-0.25, 0.0, 0.25]
    colors  = [c.color for c in cases]
    labels  = [c.label for c in cases]

    ax1 = Axis(fig_summary[1, 1],
        xlabel="Geometry", ylabel="Max δ (in)", title="Peak deflection by case",
        xticks=(1:ng, summary_df.geometry), xticklabelrotation=π/4,
        titlesize=fontsize, xlabelsize=smallfontsize, ylabelsize=smallfontsize,
        xticklabelsize=smallfontsize, yticklabelsize=smallfontsize,
    )
    δ_cols = [:max_δ_bare, :max_δ_bare25, :max_δ_composite]
    for ci in 1:3
        barplot!(ax1, collect(1:ng) .+ offsets[ci], summary_df[!, δ_cols[ci]],
            color=colors[ci], width=0.23, label=labels[ci])
    end
    axislegend(ax1, position=:rt, labelsize=smallfontsize)

    ax2 = Axis(fig_summary[2, 1],
        xlabel="Geometry", ylabel="Steel mass (kg)", title="Beam steel mass by case",
        xticks=(1:ng, summary_df.geometry), xticklabelrotation=π/4,
        titlesize=fontsize, xlabelsize=smallfontsize, ylabelsize=smallfontsize,
        xticklabelsize=smallfontsize, yticklabelsize=smallfontsize,
    )
    m_cols = [:mass_bare, :mass_bare25, :mass_composite]
    for ci in 1:3
        barplot!(ax2, collect(1:ng) .+ offsets[ci], summary_df[!, m_cols[ci]],
            color=colors[ci], width=0.23, label=labels[ci])
    end
    axislegend(ax2, position=:rt, labelsize=smallfontsize)

    save(joinpath(save_path, "composite_summary.pdf"), fig_summary)
end

# ── scatter: three panels ─────────────────────────────────────────────────────
scatter_pairs = [
    (1, 2, "Bare steel",        "Bare steel (2.5×)"),
    (1, 3, "Bare steel",        "Composite"),
    (2, 3, "Bare steel (2.5×)", "Composite"),
]

has_data = any(v -> length(v) >= 2, all_max_δ)
if has_data
    fig_scatter = Figure(size=(1200, 420), fontsize=fontsize)
    len_max = isempty(all_group_lengths) ? 1.0 : max(maximum(all_group_lengths), 1e-9)
    len_range = (0.0, len_max)
    len_cmap = :blues

    for (pi, (ci_x, ci_y, lab_x, lab_y)) in enumerate(scatter_pairs)
        δx = all_max_δ[ci_x]
        δy = all_max_δ[ci_y]
        nz = findall(i -> δx[i] > 0 || δy[i] > 0, eachindex(δx))
        length(nz) < 2 && continue

        δx_nz = δx[nz]
        δy_nz = δy[nz]
        L_nz = all_group_lengths[nz]

        ax = Axis(fig_scatter[1, pi],
            xlabel="Max δ — $lab_x (in)", ylabel="Max δ — $lab_y (in)",
            xlabelsize=smallfontsize, ylabelsize=smallfontsize,
            xticklabelsize=smallfontsize, yticklabelsize=smallfontsize,
            aspect=AxisAspect(1),
        )
        scatter!(ax, δx_nz, δy_nz; color=L_nz, colormap=len_cmap, colorrange=len_range, markersize=5, label="Groups")

        X = [ones(length(δx_nz)) δx_nz]
        β = X \ δy_nz
        lim = max(maximum(δx_nz), maximum(δy_nz)) * 1.05
        xl = range(0.0, lim; length=50)
        yl = β[1] .+ β[2] .* xl
        lines!(ax, xl, yl; color=色[:magenta], linewidth=2,
               label="Best fit (slope=$(round(β[2], digits=3)))")
        lines!(ax, [0, lim], [0, lim]; color=(:black, 0.4), linestyle=:dash,
               linewidth=1, label="1:1")
        axislegend(ax; position=:lt, labelsize=smallfontsize - 1)
    end

    Colorbar(fig_scatter[1, 4], colormap=len_cmap, colorrange=len_range,
             label="Group chain length (in)", labelsize=smallfontsize, ticklabelsize=smallfontsize - 1)
    Label(fig_scatter[0, :], "Per-beam max deflection scatter",
          fontsize=fontsize, font=:bold, halign=:center)
    save(joinpath(save_path, "composite_scatter.pdf"), fig_scatter)
end

# ── scatter: three-panel mass comparison (total + per-beam) ─────────────────────
mass_pairs = [
    (1, 2, "Bare steel",        "Bare steel (2.5×)"),
    (1, 3, "Bare steel",        "Composite"),
    (2, 3, "Bare steel (2.5×)", "Composite"),
]

if nrow(summary_df) >= 2
    L_tot = summary_df.total_beam_length
    len_max = max(maximum(L_tot), 1e-9)
    len_range = (0.0, len_max)
    len_cmap = :blues
    m_cols = [:mass_bare, :mass_bare25, :mass_composite]

    fig_mass = Figure(size=(1200, 420), fontsize=fontsize)
    for (pi, (ci_x, ci_y, lab_x, lab_y)) in enumerate(mass_pairs)
        mx = summary_df[!, m_cols[ci_x]]
        my = summary_df[!, m_cols[ci_y]]
        nz_m = findall(i -> mx[i] > 0 || my[i] > 0, eachindex(mx))
        length(nz_m) < 2 && continue

        mx_nz = mx[nz_m]
        my_nz = my[nz_m]
        L_nz = L_tot[nz_m]

        ax = Axis(fig_mass[1, pi],
            xlabel="Mass — $lab_x (kg)", ylabel="Mass — $lab_y (kg)",
            xlabelsize=smallfontsize, ylabelsize=smallfontsize,
            xticklabelsize=smallfontsize, yticklabelsize=smallfontsize,
            aspect=AxisAspect(1),
        )
        scatter!(ax, mx_nz, my_nz; color=L_nz, colormap=len_cmap, colorrange=len_range,
                 markersize=5, label="Geometries")

        X = [ones(length(mx_nz)) mx_nz]
        β = X \ my_nz
        lim = max(maximum(mx_nz), maximum(my_nz)) * 1.05
        xl = range(0.0, lim; length=50)
        yl = β[1] .+ β[2] .* xl
        lines!(ax, xl, yl; color=色[:magenta], linewidth=2,
               label="Best fit (slope=$(round(β[2], digits=3)))")
        lines!(ax, [0, lim], [0, lim]; color=(:black, 0.4), linestyle=:dash,
               linewidth=1, label="1:1")
        axislegend(ax; position=:lt, labelsize=smallfontsize - 1)
    end

    Colorbar(fig_mass[1, 4], colormap=len_cmap, colorrange=len_range,
             label="Total beam length (m)", labelsize=smallfontsize, ticklabelsize=smallfontsize - 1)
    Label(fig_mass[0, :], "Design mass comparison (total)",
          fontsize=fontsize, font=:bold, halign=:center)
    save(joinpath(save_path, "composite_mass_scatter.pdf"), fig_mass)
end

# ── scatter: per-beam mass comparison ─────────────────────────────────────────
has_mass_data = any(v -> length(v) >= 2, all_group_mass)
if has_mass_data
    fig_perbeam = Figure(size=(1200, 420), fontsize=fontsize)
    len_max = isempty(all_group_lengths) ? 1.0 : max(maximum(all_group_lengths), 1e-9)
    len_range = (0.0, len_max)
    len_cmap = :blues

    for (pi, (ci_x, ci_y, lab_x, lab_y)) in enumerate(mass_pairs)
        mx = all_group_mass[ci_x]
        my = all_group_mass[ci_y]
        nz = findall(i -> mx[i] > 0 || my[i] > 0, eachindex(mx))
        length(nz) < 2 && continue

        mx_nz = mx[nz]
        my_nz = my[nz]
        L_nz = all_group_lengths[nz]

        ax = Axis(fig_perbeam[1, pi],
            xlabel="Mass — $lab_x (kg)", ylabel="Mass — $lab_y (kg)",
            xlabelsize=smallfontsize, ylabelsize=smallfontsize,
            xticklabelsize=smallfontsize, yticklabelsize=smallfontsize,
            aspect=AxisAspect(1),
        )
        scatter!(ax, mx_nz, my_nz; color=L_nz, colormap=len_cmap, colorrange=len_range,
                 markersize=5, label="Groups")

        X = [ones(length(mx_nz)) mx_nz]
        β = X \ my_nz
        lim = max(maximum(mx_nz), maximum(my_nz)) * 1.05
        xl = range(0.0, lim; length=50)
        yl = β[1] .+ β[2] .* xl
        lines!(ax, xl, yl; color=色[:magenta], linewidth=2,
               label="Best fit (slope=$(round(β[2], digits=3)))")
        lines!(ax, [0, lim], [0, lim]; color=(:black, 0.4), linestyle=:dash,
               linewidth=1, label="1:1")
        axislegend(ax; position=:lt, labelsize=smallfontsize - 1)
    end

    Colorbar(fig_perbeam[1, 4], colormap=len_cmap, colorrange=len_range,
             label="Group chain length (in)", labelsize=smallfontsize, ticklabelsize=smallfontsize - 1)
    Label(fig_perbeam[0, :], "Per-beam mass comparison",
          fontsize=fontsize, font=:bold, halign=:center)
    save(joinpath(save_path, "composite_perbeam_mass_scatter.pdf"), fig_perbeam)
end

println("\nDone. Figures saved to $save_path")
