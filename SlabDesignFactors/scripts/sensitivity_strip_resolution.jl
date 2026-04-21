"""
Strip-resolution sensitivity study.

For every JSON layout in Geometries/topology/, runs the full
analyze → size → postprocess pipeline at several strip spacings
(the `spacing` parameter in SlabAnalysisParams, default 0.1 m).

Collects embodied-carbon metrics, beam masses, slab depths, and max
deflections into a DataFrame and saves to CSV. This directly addresses
R1 #25, R2 #8 (discretisation convergence).
"""

include("_scripts.jl")
using PrettyTables

CairoMakie.activate!()

# ── configuration ─────────────────────────────────────────────────────────────
main_path  = "Geometries/topology/"
save_path  = "SlabDesignFactors/results/sensitivity_strip_resolution/"
mkpath(save_path)

spacings = [1.0, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001]   # metres – coarse → fine

json_files = filter(x -> endswith(x, ".json"), readdir(main_path))
if isempty(json_files)
    @warn "No JSON files found in $main_path"
end

# ── results DataFrame ─────────────────────────────────────────────────────────
results_df = DataFrame(
    geometry            = String[],
    spacing             = Float64[],
    n_cells             = Int[],       # number of cells (graph faces) in layout
    n_beams             = Int[],
    n_strips            = Int[],
    mass_beams_kg       = Float64[],
    norm_mass_beams     = Float64[],   # kg/m²
    ec_steel            = Float64[],   # kgCO₂e/m²
    ec_slab             = Float64[],   # kgCO₂e/m²
    ec_rebar            = Float64[],   # kgCO₂e/m²
    ec_total            = Float64[],   # kgCO₂e/m²
    max_deflection_in   = Float64[],
    slab_area_m2        = Float64[],
    mean_slab_depth_m   = Float64[],
    elapsed_s           = Float64[],
)

# ── main loop ─────────────────────────────────────────────────────────────────
for json_file in json_files

    name = replace(json_file, ".json" => "")
    path = joinpath(main_path, json_file)

    println("\n════════════════════════════════════════════")
    println("  $name")
    println("════════════════════════════════════════════")

    raw = JSON.parse(JSON.parse(replace(read(path, String), "\\n" => ""), dicttype=Dict))
    geometry_dict = raw isa Dict ? raw : Dict(pairs(raw))

    # Number of cells (graph faces) — same for all spacings
    geom_once, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)
    slab_params_cells = SlabAnalysisParams(geom_once, slab_name=name, slab_type=:isotropic,
        vector_1d=[1.0, 0.0], slab_sizer=:uniform, spacing=0.1, plot_analysis=false,
        fix_param=true, slab_units=:m)
    slab_nodes = [node for node in slab_params_cells.model.nodes if node.id != :fixed]
    adjacency_dict, half_edges, slab_nodes = get_half_edges(slab_nodes, slab_params_cells.model.elements)
    cycles = get_cycles(adjacency_dict, half_edges)
    valid_cycles = filter_valid_cycles(slab_params_cells, cycles)
    n_cells = length(valid_cycles)

    for s in spacings

        println("  spacing = $s m")
        t0 = time()

        try
            geom, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)

            slab_params = SlabAnalysisParams(
                geom,
                slab_name       = name,
                slab_type       = :isotropic,
                vector_1d       = [1.0, 0.0],
                slab_sizer      = :uniform,
                spacing         = s,
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
                collinear              = true,
                minimum_continuous     = true,
                n_max_sections         = 0,
                composite_action       = true,
                deflection_reduction_factor = 1.0,
            )

            slab_params = analyze_slab(slab_params)
            slab_params, sizing_params = optimal_beamsizer(slab_params, sizing_params)
            results = postprocess_slab(slab_params, sizing_params, check_collinear=false)

            dt = time() - t0

            n_beams = length(results.Δ_local)
            max_δ = maximum(maximum(abs.(d)) for d in results.Δ_local if !isempty(d); init=0.0)
            n_strips = length(slab_params.load_areas)
            mean_depth = isempty(slab_params.slab_depths) ? 0.0 : mean(slab_params.slab_depths)

            push!(results_df, (
                name,
                s,
                n_cells,
                n_beams,
                n_strips,
                results.mass_beams,
                results.norm_mass_beams,
                results.embodied_carbon_beams,
                results.embodied_carbon_slab,
                results.embodied_carbon_rebar,
                results.embodied_carbon_beams + results.embodied_carbon_slab + results.embodied_carbon_rebar,
                max_δ,
                results.area,
                mean_depth,
                dt,
            ))

            println("    ✓ $(round(dt, digits=1))s  |  EC_total=$(round(results.embodied_carbon_beams + results.embodied_carbon_slab + results.embodied_carbon_rebar, digits=3))  |  max δ=$(round(max_δ, digits=4)) in  |  strips=$n_strips")

        catch e
            dt = time() - t0
            @warn "  ✗ $name @ spacing=$s failed" exception=(e, catch_backtrace())
            push!(results_df, (name, s, n_cells, 0, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, dt))
        end
    end
end

# ── save ──────────────────────────────────────────────────────────────────────
csv_path = joinpath(save_path, "strip_resolution_sensitivity.csv")
CSV.write(csv_path, results_df)
println("\n✓ Results saved to $csv_path  ($(nrow(results_df)) rows)")

# ── convergence summary table (per geometry; true finest spacing baseline) ────
println("\n──── Convergence summary (% change relative to finest spacing = $(minimum(spacings)) m) ────")
for name in unique(results_df.geometry)
    df_g = filter(row -> row.geometry == name, results_df)
    sort!(df_g, :spacing)
    finest_row = df_g[argmin(df_g.spacing), :]  # true finest spacing
    if isnan(finest_row.ec_total)
        println("  $name: finest spacing ($(minimum(spacings)) m) failed")
        continue
    end
    println("\n  $name")
    println("  " * rpad("spacing", 10) * rpad("EC_total", 12) * rpad("Δ vs fine", 12) * rpad("max δ (in)", 12) * rpad("Δ vs fine", 12) * rpad("time (s)", 12) * "strips")
    for row in eachrow(df_g)
        if isnan(row.ec_total)
            println("  " * rpad(row.spacing, 10) * "FAILED")
            continue
        end
        ec_pct = finest_row.ec_total > 0 ? (row.ec_total - finest_row.ec_total) / finest_row.ec_total * 100 : 0.0
        δ_pct  = finest_row.max_deflection_in > 0 ? (row.max_deflection_in - finest_row.max_deflection_in) / finest_row.max_deflection_in * 100 : 0.0
        println("  " *
            rpad(row.spacing, 10) *
            rpad(round(row.ec_total, digits=3), 12) *
            rpad("$(round(ec_pct, digits=2))%", 12) *
            rpad(round(row.max_deflection_in, digits=4), 12) *
            rpad("$(round(δ_pct, digits=2))%", 12) *
            rpad(round(row.elapsed_s, digits=2), 12) *
            string(row.n_strips))
    end
end

# ── convergence plot (colored by # cells; includes runtime subplot) ────────────
cell_cmap = Reverse(:dense)   # dark blue = many cells, light blue = few cells
df_plot = filter(row -> !isnan(row.ec_total), results_df)
cell_range = isempty(df_plot) ? (0, 1) : (minimum(df_plot.n_cells), maximum(df_plot.n_cells))

fig = Figure(size=(1300, 500), fontsize=11)
ax1 = Axis(fig[1, 1],
    xlabel = "Strip spacing (m)",
    ylabel = "Total EC (kgCO₂e/m²)",
    title  = "Embodied carbon vs strip resolution",
    xscale = log10,
)
ax2 = Axis(fig[1, 2],
    xlabel = "Strip spacing (m)",
    ylabel = "Max deflection (in)",
    title  = "Max deflection vs strip resolution",
    xscale = log10,
)
ax3 = Axis(fig[1, 3],
    xlabel = "Strip spacing (m)",
    ylabel = "Runtime (s)",
    title  = "Runtime vs strip resolution",
    xscale = log10,
    yscale = log10,
)

for name in unique(results_df.geometry)
    df_g = filter(row -> row.geometry == name && !isnan(row.ec_total), results_df)
    isempty(df_g) && continue
    sort!(df_g, :spacing)
    n_c = df_g.n_cells[1]
    scatterlines!(ax1, df_g.spacing, df_g.ec_total, color=n_c, colormap=cell_cmap, colorrange=cell_range, markersize=6)
    scatterlines!(ax2, df_g.spacing, df_g.max_deflection_in, color=n_c, colormap=cell_cmap, colorrange=cell_range, markersize=6)
    scatterlines!(ax3, df_g.spacing, df_g.elapsed_s, color=n_c, colormap=cell_cmap, colorrange=cell_range, markersize=6)
end

Colorbar(fig[1, 4], colormap=cell_cmap, colorrange=cell_range, label="# cells in layout", labelsize=9, ticklabelsize=8)
save(joinpath(save_path, "strip_resolution_convergence.pdf"), fig)
println("\n✓ Plot saved to $(joinpath(save_path, "strip_resolution_convergence.pdf"))")
display(fig)