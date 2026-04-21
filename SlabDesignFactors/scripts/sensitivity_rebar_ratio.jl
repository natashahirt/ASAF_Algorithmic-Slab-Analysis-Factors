"""
Deterministic post-hoc sensitivity to the reinforcement ratio.

For each layout in the `yesdeflection_yesslabmin` result set, we perturb the
rebar mass by a fractional factor δ ∈ {-0.50, -0.25, 0, +0.25, +0.50}. To
preserve total slab volume, the additional/removed rebar *volume* is swapped
against concrete (rebar displaces concrete one-for-one by volume):

    new_rebar_norm    = rebar_norm * (1 + δ)
    new_concrete_norm = concrete_norm − δ * rebar_norm * (ρ_CONCRETE / ρ_REBAR)

EC coefficients (`ECC_STEEL`, `ECC_CONCRETE`, `ECC_REBAR`) are held at their
nominal values so the effect isolated here is the reinforcement-ratio shift
alone (no MC). Use `sensitivity_monte_carlo_ec.jl` for the coefficient MC
study; the two sensitivities are reported separately in the paper to keep
each contribution interpretable.

Total EC accounting is steel + concrete + rebar only (no fireproofing or
columns), matching the underlying CSVs.

Outputs (`SlabDesignFactors/results/sensitivity_rebar_ratio/`):
  - `rebar_ratio_perturbation.csv`  — per-layout × per-δ masses and total EC
  - `rebar_ratio_ec.pdf`            — EC vs layout across δ levels
  - `rebar_ratio_bau_reduction.pdf` — best-vs-BaU reduction (%) as a function of δ

Addresses the "uniform material properties" reviewer comment on reinforcement
ratio variation.
"""

include("_scripts.jl")

CairoMakie.activate!()

# ── configuration ─────────────────────────────────────────────────────────────
deflection   = "yes"
slabmin      = "yes"
results_base = "SlabDesignFactors/results/processed_$(deflection)deflection_$(slabmin)slabmin/"
save_path    = "SlabDesignFactors/results/sensitivity_rebar_ratio/"
mkpath(save_path)

# Fractional perturbation of rebar mass. δ = +0.25 means +25% rebar.
const δ_levels = [-0.50, -0.25, 0.0, 0.25, 0.50]

# Volume-to-mass conversion: mass of concrete displaced per unit rebar mass added
# α = ρ_CONCRETE / ρ_REBAR ≈ 2400 / 7850 ≈ 0.306
const α = ρ_CONCRETE / ρ_REBAR

# ── helpers ───────────────────────────────────────────────────────────────────
"""
    perturb_rebar(df, δ) -> NamedTuple

Returns perturbed `rebar_norm`, `concrete_norm`, and `total_ec` vectors for the
given fractional rebar-mass change δ. Steel mass is unchanged. Concrete mass is
reduced by the volume equivalent of the added rebar, clipped at zero.
"""
function perturb_rebar(df, δ::Real)
    new_rebar    = df.rebar_norm    .* (1 + δ)
    new_concrete = max.(df.concrete_norm .- δ .* df.rebar_norm .* α, 0.0)
    total_ec     = df.steel_norm .* ECC_STEEL .+
                   new_concrete  .* ECC_CONCRETE .+
                   new_rebar     .* ECC_REBAR
    return (rebar_norm = new_rebar,
            concrete_norm = new_concrete,
            total_ec = total_ec)
end

"""Locate the business-as-usual layout (r1c2 uniaxial discrete collinear)."""
function find_bau_idx(df)
    return findfirst(
        row -> row.name == "r1c2" &&
               row.slab_type == "uniaxial" &&
               row.beam_sizer == "discrete" &&
               row.collinear == true,
        eachrow(df),
    )
end

# ── load data ─────────────────────────────────────────────────────────────────
df_all = assemble_data(results_base)
df     = filter(row -> row.area > 0 && !isnan(row.steel_norm), df_all)
println("\nLoaded $(nrow(df)) valid layout configurations from $(results_base)")
println("ρ_CONCRETE / ρ_REBAR = $(round(α, digits=4))")

# ── long-format perturbation table ────────────────────────────────────────────
records = DataFrame(
    name          = String[],
    rowcol        = String[],
    category      = String[],
    slab_type     = Any[],
    beam_sizer    = Any[],
    collinear     = Bool[],
    delta         = Float64[],
    rebar_norm    = Float64[],
    concrete_norm = Float64[],
    total_ec      = Float64[],
)
for δ in δ_levels
    p = perturb_rebar(df, δ)
    for i in 1:nrow(df)
        push!(records, (
            name          = df.name[i],
            rowcol        = df.rowcol[i],
            category      = df.category[i],
            slab_type     = df.slab_type[i],
            beam_sizer    = df.beam_sizer[i],
            collinear     = df.collinear[i],
            delta         = δ,
            rebar_norm    = p.rebar_norm[i],
            concrete_norm = p.concrete_norm[i],
            total_ec      = p.total_ec[i],
        ))
    end
end

csv_path = joinpath(save_path, "rebar_ratio_perturbation.csv")
CSV.write(csv_path, records)
println("✓ Per-layout × per-δ results → $csv_path  ($(nrow(records)) rows)")

# ── BaU reduction table and best-layout identity vs δ ────────────────────────
bau_idx = find_bau_idx(df)

println("\n── Reduction best-vs-BaU by rebar-ratio perturbation ──")
reductions  = Float64[]
best_rowcol = String[]
if !isnothing(bau_idx)
    for δ in δ_levels
        p         = perturb_rebar(df, δ)
        bau_ec    = p.total_ec[bau_idx]
        best_idx  = argmin(p.total_ec)
        reduction = (bau_ec - p.total_ec[best_idx]) / bau_ec * 100
        push!(reductions, reduction)
        push!(best_rowcol, df.rowcol[best_idx])
        println("  δ = $(lpad(string(round(Int, δ*100)), 4))%   BaU EC=$(round(bau_ec, digits=2))   best=$(df.rowcol[best_idx])   EC=$(round(p.total_ec[best_idx], digits=2))   reduction=$(round(reduction, digits=1))%")
    end
else
    println("  (BaU index not found — skipping reduction table)")
end

# ── plots ─────────────────────────────────────────────────────────────────────
色 = Dict(
    :ceruleanblue => colorant"#00AEEF",
    :charcoalgrey => colorant"#3e3e3e",
    :magenta      => colorant"#dc267f",
    :gold         => colorant"#df7e00",
    :lilac        => colorant"#A678B5",
)
FONTSIZE      = 11
SMALLFONTSIZE = 8

# (1) EC vs layout across δ levels, sorted by δ=0 EC
let
    nominal_sub = filter(row -> row.delta == 0.0, records)
    sort!(nominal_sub, :total_ec)
    n_show        = min(nrow(nominal_sub), 60)
    ordered_names = nominal_sub.rowcol[1:n_show]

    fig = Figure(size = (1100, 500), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "Layout (sorted by δ=0 total EC)",
        ylabel = "Total EC (kgCO₂e/m²)",
        title  = "EC sensitivity to rebar-mass perturbation (volume swap, nominal coefficients)",
        xticklabelsize = 6,
        xticklabelrotation = π/2,
        xticks = (1:n_show, ordered_names),
    )
    palette = [色[:ceruleanblue], 色[:charcoalgrey], 色[:magenta], 色[:gold], 色[:lilac]]
    for (k, δ) in enumerate(δ_levels)
        sub       = filter(row -> row.delta == δ, records)
        ec_lookup = Dict(sub.rowcol .=> sub.total_ec)
        y         = [ec_lookup[nm] for nm in ordered_names]
        scatter!(ax, 1:n_show, y,
            color      = (palette[mod1(k, length(palette))], 0.8),
            markersize = 4,
            label      = "δ = $(round(Int, δ*100))%")
    end
    axislegend(ax, position = :lt, labelsize = SMALLFONTSIZE, framevisible = false)
    save(joinpath(save_path, "rebar_ratio_ec.pdf"), fig)
    println("\n✓ EC-vs-δ plot saved")
end

# (2) Best-vs-BaU reduction as a function of δ
if !isnothing(bau_idx)
    fig = Figure(size = (600, 400), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "Rebar-mass perturbation δ (%)",
        ylabel = "EC reduction best-vs-BaU (%)",
        title  = "Best-vs-BaU reduction sensitivity to rebar ratio",
    )
    lines!(ax,    δ_levels .* 100, reductions,
        color = 色[:ceruleanblue], linewidth = 2)
    scatter!(ax, δ_levels .* 100, reductions,
        color = 色[:magenta], markersize = 8)
    save(joinpath(save_path, "rebar_ratio_bau_reduction.pdf"), fig)
    println("✓ BaU-reduction-vs-δ plot saved")
end

println("\n════════════════════════════════════════════")
println("  Rebar ratio sensitivity complete")
println("  Levels: $(round.(Int, δ_levels .* 100))%")
println("  Results: $save_path")
println("════════════════════════════════════════════")
