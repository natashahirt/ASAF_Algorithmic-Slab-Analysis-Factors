using Pkg
const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(_REPO_ROOT)

using CSV
using DataFrames
using Statistics

"""
    select_best_rows(df) -> DataFrame

For each (geometry, solver), keep one representative row:
1) prefer successful rows over failed rows
2) within that subset, keep the last row in file order
"""
function select_best_rows(df::DataFrame)::DataFrame
    required = [:geometry, :solver, :success]
    @assert all(c -> c in names(df), required) "Input is missing required columns: $(required)"

    buckets = Dict{Tuple{String, String}, Vector{Int}}()
    for (i, row) in enumerate(eachrow(df))
        key = (String(row.geometry), String(row.solver))
        get!(buckets, key, Int[])
        push!(buckets[key], i)
    end

    chosen = Int[]
    for (_, idxs) in buckets
        success_idxs = [i for i in idxs if Bool(df.success[i])]
        push!(chosen, isempty(success_idxs) ? last(idxs) : last(success_idxs))
    end
    sort!(chosen)
    return df[chosen, :]
end

"""
    safe_pct_delta(x, x_ref)

Return percentage delta `(x - x_ref) / x_ref * 100`, handling NaN and zero baseline.
"""
function safe_pct_delta(x::Real, x_ref::Real)::Float64
    if isnan(x) || isnan(x_ref) || x_ref == 0
        return NaN
    end
    return (x - x_ref) / x_ref * 100
end

"""
    safe_ratio(x, x_ref)

Return ratio `x / x_ref`, handling NaN and zero baseline.
"""
function safe_ratio(x::Real, x_ref::Real)::Float64
    if isnan(x) || isnan(x_ref) || x_ref == 0
        return NaN
    end
    return x / x_ref
end

"""
    finite_values(x)

Collect finite values from a vector-like column.
"""
function finite_values(x)
    vals = Float64[]
    for v in x
        if ismissing(v)
            continue
        end
        fv = Float64(v)
        if isfinite(fv)
            push!(vals, fv)
        end
    end
    return vals
end

safe_mean(x) = (vals = finite_values(x); isempty(vals) ? NaN : mean(vals))
safe_median(x) = (vals = finite_values(x); isempty(vals) ? NaN : median(vals))

"""
    build_comparison_tables(df_best)

Build:
1) per-geometry solver vs MIP comparison table
2) per-solver aggregated summary
"""
function build_comparison_tables(df_best::DataFrame)
    required = [:geometry, :solver, :success, :elapsed_s, :norm_mass_beams, :ec_total, :max_deflection_in]
    @assert all(c -> c in names(df_best), required) "Input is missing required columns: $(required)"

    # Baselines by geometry from successful MIP rows.
    mip_by_geometry = Dict{String, NamedTuple}()
    for row in eachrow(df_best)
        if String(row.solver) == "MIP" && Bool(row.success)
            mip_by_geometry[String(row.geometry)] = (
                elapsed_s=Float64(row.elapsed_s),
                norm_mass_beams=Float64(row.norm_mass_beams),
                ec_total=Float64(row.ec_total),
                max_deflection_in=Float64(row.max_deflection_in),
            )
        end
    end

    comparison_df = DataFrame(
        geometry=String[],
        solver=String[],
        success=Bool[],
        elapsed_s=Float64[],
        norm_mass_beams=Float64[],
        ec_total=Float64[],
        max_deflection_in=Float64[],
        runtime_ratio_vs_mip=Float64[],
        norm_mass_delta_pct_vs_mip=Float64[],
        ec_total_delta_pct_vs_mip=Float64[],
        max_deflection_delta_pct_vs_mip=Float64[],
    )

    for row in eachrow(df_best)
        geom = String(row.geometry)
        solver = String(row.solver)
        base = get(mip_by_geometry, geom, nothing)

        elapsed_s = Float64(row.elapsed_s)
        mass = Float64(row.norm_mass_beams)
        ec = Float64(row.ec_total)
        defl = Float64(row.max_deflection_in)

        runtime_ratio = isnothing(base) ? NaN : safe_ratio(elapsed_s, base.elapsed_s)
        mass_delta = isnothing(base) ? NaN : safe_pct_delta(mass, base.norm_mass_beams)
        ec_delta = isnothing(base) ? NaN : safe_pct_delta(ec, base.ec_total)
        defl_delta = isnothing(base) ? NaN : safe_pct_delta(defl, base.max_deflection_in)

        push!(comparison_df, (
            geom,
            solver,
            Bool(row.success),
            elapsed_s,
            mass,
            ec,
            defl,
            runtime_ratio,
            mass_delta,
            ec_delta,
            defl_delta,
        ))
    end

    # Aggregate by solver, excluding MIP (baseline).
    non_mip = filter(r -> r.solver != "MIP", comparison_df)
    solver_summary = combine(groupby(non_mip, :solver),
        :success => (x -> sum(x)) => :n_success,
        :success => length => :n_total,
        :runtime_ratio_vs_mip => safe_mean => :mean_runtime_ratio_vs_mip,
        :runtime_ratio_vs_mip => safe_median => :median_runtime_ratio_vs_mip,
        :ec_total_delta_pct_vs_mip => safe_mean => :mean_ec_delta_pct_vs_mip,
        :ec_total_delta_pct_vs_mip => safe_median => :median_ec_delta_pct_vs_mip,
        :norm_mass_delta_pct_vs_mip => safe_mean => :mean_mass_delta_pct_vs_mip,
        :max_deflection_delta_pct_vs_mip => safe_mean => :mean_deflection_delta_pct_vs_mip,
    )
    solver_summary.success_rate_pct = 100 .* solver_summary.n_success ./ solver_summary.n_total

    return comparison_df, solver_summary
end

function main()
    args = ARGS
    if !(length(args) in (1, 2))
        println("Usage: julia executable_nlp_solver_summary.jl <results_path> [output_path]")
        println("Reads <results_path>/nlp_solver_comparison.csv by default.")
        return
    end

    results_path = args[1]
    output_path = length(args) == 2 ? args[2] : results_path

    input_csv = joinpath(results_path, "nlp_solver_comparison.csv")
    if !isfile(input_csv)
        error("Input file not found: $input_csv")
    end

    mkpath(output_path)

    raw_df = CSV.read(input_csv, DataFrame)
    best_df = select_best_rows(raw_df)
    comparison_df, solver_summary = build_comparison_tables(best_df)

    per_geometry_csv = joinpath(output_path, "nlp_solver_vs_mip_per_geometry.csv")
    summary_csv = joinpath(output_path, "nlp_solver_vs_mip_summary.csv")
    CSV.write(per_geometry_csv, comparison_df)
    CSV.write(summary_csv, solver_summary)

    println("Wrote per-geometry comparison: $per_geometry_csv")
    println("Wrote solver summary: $summary_csv")
end

main()
