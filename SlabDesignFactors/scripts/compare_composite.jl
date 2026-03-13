"""
Compare beam deflections with and without composite action for all topology geometries.

Runs each JSON in Geometries/topology/ through:
  1. Bare steel (composite_action=false, the existing default)
  2. Composite action (composite_action=true, AISC transformed section)

Produces per-geometry figures overlaying the deflection envelope for every beam,
and a summary table of max deflections and steel weight for both cases.
"""

include("_scripts.jl")

CairoMakie.activate!()

# ── colour palette (matches existing plot files) ──────────────────────────────
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

color_bare      = 色[:charcoalgrey]
color_composite = 色[:ceruleanblue]

fontsize      = 11
smallfontsize = 8

# ── paths ─────────────────────────────────────────────────────────────────────
main_path    = "Geometries/topology/"
save_path    = "SlabDesignFactors/plot/figures/composite/"
mkpath(save_path)

json_files = filter(x -> endswith(x, ".json"), readdir(main_path))

if isempty(json_files)
    @warn "No JSON files found in $main_path"
end

# ── summary collection ────────────────────────────────────────────────────────
summary_df = DataFrame(
    geometry       = String[],
    n_beams        = Int[],
    max_δ_bare     = Float64[],     # in
    max_δ_composite = Float64[],    # in
    mass_bare      = Float64[],     # kg
    mass_composite = Float64[],     # kg
    mass_savings_pct = Float64[],   # %
)

# ── main loop ─────────────────────────────────────────────────────────────────
for json_file in json_files

    name = replace(json_file, ".json" => "")
    path = joinpath(main_path, json_file)

    println("\n════════════════════════════════════════════")
    println("  $name")
    println("════════════════════════════════════════════")

    geometry_dict = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))

    # ── helper: run one configuration ─────────────────────────────────────
    function run_case(; composite::Bool)
        geom, type_info = generate_from_json(geometry_dict, plot=false, drawn=false)

        slab_params = SlabAnalysisParams(
            geom,
            slab_name       = name,
            slab_type       = :uniaxial,
            vector_1d       = [1.0, 0.0],
            slab_sizer      = :uniform,
            spacing         = 0.1,
            plot_analysis   = false,
            fix_param       = true,
            slab_units      = :m,
        )

        sizing_params = SlabSizingParams(
            live_load              = psf_to_ksi(50),
            superimposed_dead_load = psf_to_ksi(15),
            slab_dead_load         = 0.0,
            live_factor            = 1.6,
            dead_factor            = 1.2,
            beam_sizer             = :discrete,
            max_depth              = 40.0,
            beam_units             = :in,
            serviceability_lim     = 360,
            collinear              = false,
            minimum_continuous     = true,
            n_max_sections         = 0,
            composite_action       = composite,
        )

        slab_params = analyze_slab(slab_params)

        t0 = time()
        slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)
        dt = time() - t0
        tag = composite ? "composite" : "bare steel"
        println("  [$tag] optimal_beamsizer: $(round(dt, digits=2))s")

        results = postprocess_slab(slab_params, sizing_params, check_collinear=false)
        return results
    end

    # ── run both cases ────────────────────────────────────────────────────
    results_bare      = run_case(composite=false)
    results_composite = run_case(composite=true)

    n_beams = length(results_bare.Δ_local)

    # ── per-beam deflection comparison figure ─────────────────────────────
    n_cols = min(n_beams, 4)
    n_rows = ceil(Int, n_beams / n_cols)

    fig = Figure(size=(300 * n_cols, 220 * n_rows), fontsize=fontsize)

    max_δ_bare_all      = 0.0
    max_δ_composite_all = 0.0

    for i in 1:n_beams
        row = div(i - 1, n_cols) + 1
        col = mod(i - 1, n_cols) + 1

        ax = Axis(fig[row, col],
            xlabel = "x (in)", ylabel = "δ (in)",
            title  = "Beam $i: $(results_bare.ids[i])",
            titlesize = smallfontsize,
            xlabelsize = smallfontsize, ylabelsize = smallfontsize,
            xticklabelsize = smallfontsize, yticklabelsize = smallfontsize,
        )

        x_bare = results_bare.x[i]
        δ_bare = results_bare.Δ_local[i]
        x_comp = results_composite.x[i]
        δ_comp = results_composite.Δ_local[i]

        lines!(ax, x_bare, δ_bare, color=color_bare, linewidth=1.5, label="Bare steel")
        lines!(ax, x_comp, δ_comp, color=color_composite, linewidth=1.5, label="Composite")

        # serviceability limit line
        L = x_bare[end]
        δ_limit = L / 360
        hlines!(ax, [-δ_limit], color=色[:magenta], linewidth=0.8, linestyle=:dash)

        max_δ_bare_all      = max(max_δ_bare_all,      maximum(abs.(δ_bare)))
        max_δ_composite_all = max(max_δ_composite_all, maximum(abs.(δ_comp)))
    end

    # shared legend
    elem_bare = LineElement(color=color_bare, linewidth=1.5)
    elem_comp = LineElement(color=color_composite, linewidth=1.5)
    elem_lim  = LineElement(color=色[:magenta], linewidth=0.8, linestyle=:dash)
    Legend(fig[0, 1:n_cols],
        [elem_bare, elem_comp, elem_lim],
        ["Bare steel", "Composite", "L/360"],
        orientation=:horizontal, framevisible=false, labelsize=smallfontsize,
    )

    Label(fig[-1, 1:n_cols], "$name — deflection comparison",
        fontsize=fontsize, font=:bold, halign=:center)

    save(joinpath(save_path, "$(name)_deflection_comparison.pdf"), fig)
    display(fig)

    # ── summary stats ─────────────────────────────────────────────────────
    mass_savings = results_bare.mass_beams > 0 ?
        (results_bare.mass_beams - results_composite.mass_beams) / results_bare.mass_beams * 100 : 0.0

    push!(summary_df, (
        name,
        n_beams,
        max_δ_bare_all,
        max_δ_composite_all,
        results_bare.mass_beams,
        results_composite.mass_beams,
        mass_savings,
    ))

    println("  max δ bare:      $(round(max_δ_bare_all, digits=4)) in")
    println("  max δ composite: $(round(max_δ_composite_all, digits=4)) in")
    println("  steel mass bare:      $(round(results_bare.mass_beams, digits=1)) kg")
    println("  steel mass composite: $(round(results_composite.mass_beams, digits=1)) kg")
    println("  mass savings: $(round(mass_savings, digits=1))%")

    GC.gc()
end

# ── print summary table ──────────────────────────────────────────────────────
println("\n\n╔══════════════════════════════════════════════════════════════╗")
println("║           COMPOSITE ACTION COMPARISON SUMMARY              ║")
println("╚══════════════════════════════════════════════════════════════╝\n")

pretty_table(summary_df,
    header = ["Geometry", "# Beams", "max δ bare (in)", "max δ composite (in)",
              "Mass bare (kg)", "Mass composite (kg)", "Savings (%)"],
    formatters = (ft_printf("%.4f", [3, 4]), ft_printf("%.1f", [5, 6, 7])),
    crop = :none,
)

# ── summary bar chart ────────────────────────────────────────────────────────
if nrow(summary_df) > 0
    fig_summary = Figure(size=(max(600, 80 * nrow(summary_df)), 500), fontsize=fontsize)

    ax1 = Axis(fig_summary[1, 1],
        xlabel = "Geometry", ylabel = "Max δ (in)",
        title  = "Peak deflection: bare steel vs composite",
        xticks = (1:nrow(summary_df), summary_df.geometry),
        xticklabelrotation = π/4,
        titlesize = fontsize, xlabelsize = smallfontsize, ylabelsize = smallfontsize,
        xticklabelsize = smallfontsize, yticklabelsize = smallfontsize,
    )

    barplot!(ax1, collect(1:nrow(summary_df)) .- 0.15, summary_df.max_δ_bare,
        color=color_bare, width=0.3, label="Bare steel")
    barplot!(ax1, collect(1:nrow(summary_df)) .+ 0.15, summary_df.max_δ_composite,
        color=color_composite, width=0.3, label="Composite")
    axislegend(ax1, position=:rt, framevisible=true, backgroundcolor=:white,
        framecolor=:white, labelsize=smallfontsize)

    ax2 = Axis(fig_summary[2, 1],
        xlabel = "Geometry", ylabel = "Steel mass (kg)",
        title  = "Beam steel mass: bare steel vs composite",
        xticks = (1:nrow(summary_df), summary_df.geometry),
        xticklabelrotation = π/4,
        titlesize = fontsize, xlabelsize = smallfontsize, ylabelsize = smallfontsize,
        xticklabelsize = smallfontsize, yticklabelsize = smallfontsize,
    )

    barplot!(ax2, collect(1:nrow(summary_df)) .- 0.15, summary_df.mass_bare,
        color=color_bare, width=0.3, label="Bare steel")
    barplot!(ax2, collect(1:nrow(summary_df)) .+ 0.15, summary_df.mass_composite,
        color=color_composite, width=0.3, label="Composite")
    axislegend(ax2, position=:rt, framevisible=true, backgroundcolor=:white,
        framecolor=:white, labelsize=smallfontsize)

    save(joinpath(save_path, "composite_summary.pdf"), fig_summary)
    display(fig_summary)
end

println("\nDone. Figures saved to $save_path")
