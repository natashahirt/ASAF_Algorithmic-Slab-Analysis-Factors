"""
    print_forces(results::SlabOptimResults; i::Int64=-1)

Prints the forces and related data for a given slab optimization result.

# Arguments
- `results::SlabOptimResults`: The results of the slab optimization.
- `i::Int64=-1`: Index of the result to print. If -1, prints all results.

# Output
Prints a DataFrame of forces and a summary of mass and embodied carbon.
"""
function print_forces(results::SlabOptimResults; i::Int64=-1)
    df = create_forces_dataframe(results, i)
    print_forces_dataframe(df)
    print_mass_and_carbon_summary(results)
end

function create_forces_dataframe(results::SlabOptimResults, i::Int64)
    df = DataFrame(id=String[], x=Float64[], P=Float64[], My=Float64[], Mn=Float64[], Vy=Float64[], Vn=Float64[], Δ_local=Float64[], Δ_global=Float64[])
    
    for i in (i > 0 ? i : 1:lastindex(results.minimums))
        id = results.ids[i]
        x = results.x[i][end]

        # Determine maximum load, moment, shear, and deflection
        P = isempty(results.P[i]) ? 0 : maximum(abs.(results.P[i]))
        My = maximum(abs.(results.My[i]))
        Mn = results.Mn[i]

        Vy = maximum(abs.(results.Vy[i]))
        Vn = results.Vn[i]

        Δ_local = maximum(abs.(results.Δ_local[i]))
        Δ_global = maximum(abs.(results.Δ_global[i]))

        push!(df, [id, x, P, My, Mn, Vy, Vn, Δ_local, Δ_global])
    end

    return df
end

function print_forces_dataframe(df::DataFrame)
    println("Forces and Related Data:")
    println(df)
    println("\n")
end

function print_mass_and_carbon_summary(results::SlabOptimResults)
    println("Summary of Mass and Embodied Carbon:")
    println("Beams   || Steel    || Total mass: $(round(results.mass_beams, digits=2)) kg")
    println("        ||          || Normalized mass: $(round(results.norm_mass_beams, digits=2)) kg/m²")
    println("        ||          || Embodied carbon: $(round(results.embodied_carbon_beams, digits=2)) kgCO₂eq/m²")
    println("Columns || Steel    || Total mass: $(round(results.mass_columns, digits=2)) kg")
    println("        ||          || Normalized mass: $(round(results.norm_mass_columns, digits=2)) kg/m²")
    println("        ||          || Embodied carbon: $(round(results.embodied_carbon_columns, digits=2)) kgCO₂eq/m²")
    println("Slab    || Concrete || Total mass: $(round(results.mass_slab, digits=2)) kg")
    println("        ||          || Normalized mass: $(round(results.norm_mass_slab, digits=2)) kg/m²")
    println("        ||          || Embodied carbon: $(round(results.embodied_carbon_slab, digits=2)) kgCO₂eq/m²")
    println("Slab    || Rebar    || Total mass: $(round(results.mass_rebar, digits=2)) kg")
    println("        ||          || Normalized mass: $(round(results.norm_mass_rebar, digits=2)) kg/m²")
    println("        ||          || Embodied carbon: $(round(results.embodied_carbon_rebar, digits=2)) kgCO₂eq/m²")
    println("Fire-   || proofing || Total mass: $(round(results.mass_fireproofing, digits=2)) kg")
    println("        ||          || Normalized mass: $(round(results.norm_mass_fireproofing, digits=2)) kg/m²")
    println("        ||          || Embodied carbon: $(round(results.embodied_carbon_fireproofing, digits=2)) kgCO₂eq/m²")
    total_mass = results.mass_beams + results.mass_columns + results.mass_slab + results.mass_rebar + results.mass_fireproofing
    total_norm = results.norm_mass_beams + results.norm_mass_columns + results.norm_mass_slab + results.norm_mass_rebar + results.norm_mass_fireproofing
    total_ec   = results.embodied_carbon_beams + results.embodied_carbon_columns + results.embodied_carbon_slab + results.embodied_carbon_rebar + results.embodied_carbon_fireproofing
    println("Total   ||          || Total mass: $(round(total_mass, digits=2)) kg")
    println("        ||          || Normalized mass: $(round(total_norm, digits=2)) kg/m²")
    println("        ||          || Embodied carbon: $(round(total_ec, digits=2)) kgCO₂eq/m²")
    if !isempty(results.col_sections)
        println("\nColumn Sizing:")
        for i in 1:length(results.col_sections)
            println("  Col $i: $(results.col_sections[i])  " *
                    "Pu=$(round(results.col_Pu[i], digits=1)) kip  " *
                    "ϕPn=$(round(results.col_ϕPn[i], digits=1)) kip  " *
                    "util=$(round(results.col_util[i]*100, digits=1))%")
        end
    end
    if !isempty(results.δ_total)
        has_beam_dead = !isempty(results.δ_beam_dead)
        println("\nStaged Deflection (unshored composite):")
        hdr = "  Beam  │  δ_DL(bare)  │"
        hdr *= has_beam_dead ? "  δ_beam_SW  │" : ""
        hdr *= "  δ_SDL(comp)  │  δ_LL(comp)  │  δ_total  │  L/360  │  L/240  │ LL ok │ Tot ok"
        println(hdr)
        IN_TO_MM_TBL = 25.4
        for i in 1:length(results.δ_total)
            id = i <= length(results.ids) ? results.ids[i] : "$i"
            row = "  $(lpad(id, 5)) │ " *
                  "$(lpad(round(results.δ_slab_dead[i]*IN_TO_MM_TBL, digits=2), 11)) │ "
            if has_beam_dead
                row *= "$(lpad(round(results.δ_beam_dead[i]*IN_TO_MM_TBL, digits=2), 10)) │ "
            end
            row *= "$(lpad(round(results.δ_sdl[i]*IN_TO_MM_TBL, digits=2), 12)) │ " *
                   "$(lpad(round(results.δ_live[i]*IN_TO_MM_TBL, digits=2), 11)) │ " *
                   "$(lpad(round(results.δ_total[i]*IN_TO_MM_TBL, digits=2), 8)) │ " *
                   "$(lpad(round(results.Δ_limit_live[i]*IN_TO_MM_TBL, digits=2), 6)) │ " *
                   "$(lpad(round(results.Δ_limit_total[i]*IN_TO_MM_TBL, digits=2), 6)) │ " *
                   "$(lpad(results.δ_live_ok[i] ? "✓" : "✗", 5)) │ " *
                   "$(results.δ_total_ok[i] ? "✓" : "✗")"
            println(row)
        end
        println("  (all deflections in mm)")
    end
    println("\n")
end

"""
    save_results(results_list::Vector{SlabOptimResults}; verbose::Bool=true, subfolder::String="", filename::String="")

Saves the results of slab optimizations to a CSV file.

# Arguments
- `results_list::Vector{SlabOptimResults}`: A list of slab optimization results.
- `verbose::Bool=true`: If true, prints detailed information about each result.
- `subfolder::String=""`: Subfolder to save the results in.
- `filename::String=""`: Name of the file to save the results. If empty, uses a default naming scheme.

# Output
Saves a CSV file with the results.
"""
function save_results(results_list::Vector{SlabOptimResults}; verbose::Bool=true, subfolder::String="", filename::String="")
    # Create subfolder if it doesn't exist
    if !isdir(subfolder)
        mkpath(subfolder)
    end
    df = create_results_dataframe(results_list, verbose)

    save_dataframe_to_csv(df, results_list, subfolder, filename)
end

function determine_collinear(results_list::Vector{SlabOptimResults})
    return Union{Bool,Nothing}[results_list[i].collinear for i in 1:lastindex(results_list)]
end

"""
    sync_csv_schema!(existing_df::DataFrame, new_df::DataFrame)

Ensure `existing_df` and `new_df` share the same column names so appends and row
updates work when `create_results_dataframe` gains columns (e.g. staged deflection
metadata). New columns are filled with `missing` on the side that lacked them.
"""
function sync_csv_schema!(existing_df::DataFrame, new_df::DataFrame)
    for col in names(new_df)
        if !(col in names(existing_df))
            existing_df[!, col] = Union{Missing, Any}[missing for _ in 1:nrow(existing_df)]
        end
    end
    for col in names(existing_df)
        if !(col in names(new_df))
            new_df[!, col] = Union{Missing, Any}[missing for _ in 1:nrow(new_df)]
        end
    end
    nothing
end

"""
    create_results_dataframe(results_list, verbose)

Build one row per `SlabOptimResults` for CSV export. Includes per-beam staged
deflections (mm), L/360/L/240 pass flags, failing beam indices, global deflection
sanity, strength utilization, and `optimal_beamsizer` staged-loop status — same
information as `print_mass_and_carbon_summary` / the per-beam table. The canonical
pipeline and assertions live in `SlabDesignFactors/test/run.jl` → `test/runtests.jl`
(integration tests on real topologies).
"""
function create_results_dataframe(results_list::Vector{SlabOptimResults}, verbose::Bool)
    df = DataFrame(
        name=String[], 
        area=Float64[], 
        steel_norm=Float64[],
        column_norm=Float64[],
        concrete_norm=Float64[], 
        rebar_norm=Float64[], 
        fireproofing_norm=Float64[],
        max_depth=Float64[], 
        slab_type=String[], 
        slab_sizer=String[], 
        beam_sizer=String[], 
        nlp_solver=String[],
        deflection_limit=Bool[],
        collinear=Bool[], 
        vector_1d_x=Float64[], 
        vector_1d_y=Float64[], 
        sections=String[], 
        ids=String[],
        col_sections=String[],
        col_Pu=String[],
        col_ϕPn=String[],
        col_util=String[],
        δ_slab_dead_mm=String[],
        δ_beam_dead_mm=String[],
        δ_sdl_mm=String[],
        δ_live_mm=String[],
        δ_total_mm=String[],
        δ_live_ok=String[],
        δ_total_ok=String[],
        Δ_limit_live_mm=String[],
        Δ_limit_total_mm=String[],
        composite_action=Bool[],
        staged_converged=Bool[],
        staged_n_violations=Int[],
        n_L360_fail=Int[],
        n_L240_fail=Int[],
        i_L360_fail=String[],
        i_L240_fail=String[],
        max_δ_total_mm=Float64[],
        max_bay_span_in=Float64[],
        global_δ_ok=Bool[],
        max_util_M=Float64[],
        max_util_V=Float64[],
        max_col_util=Float64[],
    )

    empty_arr = "Any[]"
    _ints(v::Vector{Int}) = isempty(v) ? "Any[]" : "Any[" * join(map(string, v), ", ") * "]"

    for i in 1:lastindex(results_list)
        results = results_list[i]
        name = results.slab_name
        area = results.area

        if area == 0
            push!(df, [name, 0., 0., 0., 0., 0., 0., results.max_depth,
                       String(results.slab_type), String(results.slab_sizer), String(results.beam_sizer),
                       results.nlp_solver, results.deflection_limit,
                       results.collinear, results.vector_1d[1], results.vector_1d[2],
                       empty_arr, empty_arr, empty_arr, empty_arr, empty_arr, empty_arr,
                       empty_arr, empty_arr, empty_arr, empty_arr, empty_arr, empty_arr, empty_arr,
                       empty_arr, empty_arr,
                       false, true, 0, 0, 0, empty_arr, empty_arr,
                       0.0, 0.0, true, 0.0, 0.0, 0.0])
            continue
        end

        string_sections = "Any[" * join(map(x -> "\"$x\"", results.sections), ", ") * "]"
        string_ids = "Any[" * join(map(x -> "\"$x\"", results.ids), ", ") * "]"
        string_col_sections = "Any[" * join(map(x -> "\"$x\"", results.col_sections), ", ") * "]"
        string_col_Pu   = "Any[" * join(map(x -> string(round(x, digits=2)), results.col_Pu), ", ") * "]"
        string_col_ϕPn  = "Any[" * join(map(x -> string(round(x, digits=2)), results.col_ϕPn), ", ") * "]"
        string_col_util = "Any[" * join(map(x -> string(round(x, digits=4)), results.col_util), ", ") * "]"

        _mm(v) = "Any[" * join(map(x -> string(round(x * 25.4, digits=3)), v), ", ") * "]"
        _bool(v) = "Any[" * join(map(x -> string(x), v), ", ") * "]"

        max_δ_tot_mm = isempty(results.δ_total) ? 0.0 : round(maximum(results.δ_total) * 25.4, digits=3)

        push!(df, [name, area, results.norm_mass_beams, results.norm_mass_columns,
                   results.norm_mass_slab, results.norm_mass_rebar, results.norm_mass_fireproofing,
                   results.max_depth,
                   String(results.slab_type), String(results.slab_sizer), String(results.beam_sizer),
                   results.nlp_solver, results.deflection_limit,
                   results.collinear, results.vector_1d[1], results.vector_1d[2],
                   string_sections, string_ids,
                   string_col_sections, string_col_Pu, string_col_ϕPn, string_col_util,
                   _mm(results.δ_slab_dead), _mm(results.δ_beam_dead),
                   _mm(results.δ_sdl), _mm(results.δ_live), _mm(results.δ_total),
                   _bool(results.δ_live_ok), _bool(results.δ_total_ok),
                   _mm(results.Δ_limit_live), _mm(results.Δ_limit_total),
                   results.composite_action, results.staged_converged, results.staged_n_violations,
                   results.n_L360_fail, results.n_L240_fail,
                   _ints(results.i_L360_fail), _ints(results.i_L240_fail),
                   max_δ_tot_mm, results.max_bay_span, results.global_δ_ok,
                   results.max_util_M, results.max_util_V, results.max_col_util])

        if verbose
            print_verbose_results(results)
        end
    end

    if !verbose
        println("Results DataFrame:")
        println(df)
    end

    return df
end

function print_verbose_results(results::SlabOptimResults)
    println("Slab $(results.slab_name):")
    println("  Type: $(results.slab_type) slab")
    println("  Vector: $(results.vector_1d)")
    println("  Slab sizer: $(results.slab_sizer)")
    println("  Beam sizer: $(results.beam_sizer)")
    println("  NLP solver / MIP: $(results.nlp_solver)  deflection_limit=$(results.deflection_limit)")
    println("  collinear: $(results.collinear)")

    if results.area == 0.
        println("  Slab span too large.\n")
    else
        println("  Beam steel normalized mass: $(round(results.norm_mass_beams, digits=2)) kg/m²")
        println("  Column steel normalized mass: $(round(results.norm_mass_columns, digits=2)) kg/m²")
        println("  Total normalized mass: $(round(results.norm_mass_beams + results.norm_mass_columns + results.norm_mass_slab + results.norm_mass_rebar, digits=2)) kg/m²")
        if !isempty(results.δ_total)
            n_live_fail  = count(.!results.δ_live_ok)
            n_total_fail = count(.!results.δ_total_ok)
            IN_TO_MM = 25.4
            max_δ_LL  = round(maximum(results.δ_live) * IN_TO_MM, digits=2)
            max_δ_tot = round(maximum(results.δ_total) * IN_TO_MM, digits=2)
            defl_str = "  Deflection: max δ_LL=$(max_δ_LL) mm, max δ_total=$(max_δ_tot) mm"
            if !isempty(results.δ_beam_dead)
                max_δ_bm = round(maximum(results.δ_beam_dead) * IN_TO_MM, digits=2)
                defl_str *= ", max δ_beam_SW=$(max_δ_bm) mm"
            end
            defl_str *= "  [L/360 fail: $(n_live_fail), L/240 fail: $(n_total_fail)]"
            println(defl_str)
        end
        if results.composite_action
            println("  Staged sizer: converged=$(results.staged_converged), " *
                    "n_violations=$(results.staged_n_violations)")
        end
        println()
    end
end

function save_dataframe_to_csv(df::DataFrame, results_list::Vector{SlabOptimResults}, subfolder::String, filename::String)
    if filename == ""
        return
    elseif filename == "default"
        filename = "$(results_list[1].slab_type)_$(results_list[1].slab_sizer)_$(results_list[1].vector_1d)_$(results_list[1].max_depth)"
    end

    folder = "results/"
    folder = ""

    if subfolder == "default"
        subfolder = "$(results_list[1].slab_type)_$(results_list[1].slab_sizer)_$(results_list[1].vector_1d)_$(results_list[1].beam_sizer)_$(results_list[1].max_depth)_$(results_list[1].collinear)"
    end

    """if occursin(folder, subfolder)
        subfolder = subfolder[length(folder)+1:end]
    end"""

    if subfolder[end] != '/' || folder == ""
        folder = folder * subfolder * "/"
    end

    filename = folder * filename * ".csv"

    println("Saving results to $folder")

    if !isdir(folder)
        mkdir(folder)
    end

    CSV.write(filename, df)
end

function append_results_to_csv(folder::String, filename::String, results_list::Vector{SlabOptimResults}; unique_sections::Int=0)
    # Create directory if it doesn't exist
    if !isdir(folder)
        mkdir(folder)
    end
    
    path = folder * filename * ".csv"
    new_df = create_results_dataframe(results_list, true)

    # Add a column for unique_sections
    new_df.unique_sections = fill(unique_sections, nrow(new_df))

    try
        if !isfile(path)
            # Create new file if it doesn't exist (write to .tmp then rename for safety on timeout)
            tmp_path = path * ".tmp"
            CSV.write(tmp_path, new_df)
            Base.Filesystem.rename(tmp_path, path)
        else
            # Read existing file
            existing_df = CSV.read(path, DataFrame)

            # Check if dataframe has any columns
            if isempty(names(existing_df))
                save_results(results_list, subfolder=folder, filename=filename)
                return
            end

            # Ensure string columns can handle longer strings
            for col in names(existing_df)
                if eltype(existing_df[!, col]) <: AbstractString
                    existing_df[!, col] = convert(Vector{String}, existing_df[!, col])
                end
            end

            sync_csv_schema!(existing_df, new_df)

            # For each row in new results
            for i in 1:nrow(new_df)
                # Check if matching row exists
                matching_rows = findall(
                    (String.(existing_df.name) .== String(new_df.name[i])) .&
                    (String.(existing_df.slab_type) .== String(new_df.slab_type[i])) .&
                    (String.(existing_df.slab_sizer) .== String(new_df.slab_sizer[i])) .&
                    (String.(existing_df.beam_sizer) .== String(new_df.beam_sizer[i])) .&
                    (existing_df.collinear .== new_df.collinear[i]) .&
                    (existing_df.vector_1d_x .== new_df.vector_1d_x[i]) .&
                    (existing_df.vector_1d_y .== new_df.vector_1d_y[i]) .&
                    (existing_df.max_depth .== new_df.max_depth[i]) .&
                    (existing_df.unique_sections .== new_df.unique_sections[i])
                )
                if length(matching_rows) == 1
                    # Convert row-by-row to ensure type compatibility
                    for col in names(existing_df)
                        existing_df[matching_rows[1], col] = new_df[i, col]
                    end
                    println("Replaced row $(matching_rows[1])")
                else
                    # Append new row
                    push!(existing_df, new_df[i, :], promote=true)
                    println("Appended row $(i)")
                end
            end

            # Write to .tmp then rename so a mid-write timeout doesn't corrupt the main file
            tmp_path = path * ".tmp"
            CSV.write(tmp_path, existing_df)
            Base.Filesystem.rename(tmp_path, path)
        end
    catch e
        error("Failed to write results to CSV file: $e")
    end

    # Verify file was written
    if !isfile(path)
        error("Failed to write results to $path")
    end
end