"""
Variance decomposition (ANOVA / η², ω²) of total embodied carbon across
design decisions, run separately on the two reviewer conditions:

    1. `yesdeflection_yesslabmin` — deflection checks on, minimum slab depth on
    2. `yesdeflection_noslabmin`  — deflection checks on, minimum slab depth off

For each condition we fit two fixed-effects linear models:

  A) Original (flat layout nuisance, preserves prior headline result)
    total_ec ~ name
             + slab_type + slab_sizer + beam_sizer + collinear + max_depth
             + all pairwise interactions between the five design decisions

  B) Supplemental (nested layout nuisance decomposition)
    total_ec ~ category
             + category:row_idx + category:col_idx + category:row_idx:col_idx
             + slab_type + slab_sizer + beam_sizer + collinear + max_depth
             + all pairwise interactions between the five design decisions

and report, for every term, the sum of squares, degrees of freedom,
F-statistic, p-value, partial η² (SS_term / SS_total) and ω² (unbiased
effect size). Model A preserves the prior "layout explains ~53%" finding.
Model B decomposes that layout component as
`category → (category×row) + (category×col) + (category×row×col)` so
between-layout variation is explicitly nested under topology family and row/
column structure. The nested decomposition is additive to (not a replacement
for) the original flat-layout result.

Because the five design factors cross evenly within each layout the within-
layout design is (near-)balanced, so Type I / II / III SS agree. The
reduced-vs-full SS implementation used here is Type II by construction
(each term tested after all other terms at the same or lower order).

A companion CSV (`master_design_matrix.csv`) is also exported per condition
so the same data can be fed into a tree-based model + SHAP analysis in
Python (`shap_analysis.py`, dropped alongside the CSV).

Outputs (per condition, in `SlabDesignFactors/results/sensitivity_variance_decomposition/<tag>/`):
  - `anova_eta_squared.csv`    — term, df, SS, MS, F, p, η², ω²
  - `anova_eta_squared.pdf`    — horizontal η² bar plot
  - `master_design_matrix.csv` — cleaned design + response for SHAP
  - `fit_summary.txt`          — sample size, R², residual DoF, etc.

Cross-condition output (top-level save_path):
  - `eta_squared_comparison.csv` — side-by-side η² table
  - `eta_squared_comparison.pdf` — grouped bar plot comparing the two conditions
"""

include("_scripts.jl")

# ── lightweight dep bootstrap ─────────────────────────────────────────────────
# GLM / StatsModels / CategoricalArrays / Distributions are not part of the
# base project. Add them into the active environment on first run so the
# script is self-contained.
using Pkg
let _declared = keys(Pkg.project().dependencies)
    for pkg in ("GLM", "StatsModels", "CategoricalArrays", "Distributions")
        pkg in _declared || Pkg.add(pkg)
    end
end

using GLM, StatsModels, CategoricalArrays, Distributions
using LinearAlgebra: I  # silence unused-warning if any downstream uses it

CairoMakie.activate!()

# ── configuration ─────────────────────────────────────────────────────────────
const CONDITIONS = [
    (tag = "yesdeflection_yesslabmin",
     label = "deflection on · slab-min on",
     base  = "SlabDesignFactors/results/processed_yesdeflection_yesslabmin/"),
    (tag = "yesdeflection_noslabmin",
     label = "deflection on · slab-min off",
     base  = "SlabDesignFactors/results/processed_yesdeflection_noslabmin/"),
]

const SAVE_ROOT = "SlabDesignFactors/results/sensitivity_variance_decomposition/"
mkpath(SAVE_ROOT)

# Design decisions whose main effects + pairwise interactions are decomposed.
const DECISIONS = [:slab_type, :slab_sizer, :beam_sizer, :collinear, :max_depth]
const LAYOUT = :name

const FONTSIZE      = 14
const SMALLFONTSIZE = 11

# ── helpers ───────────────────────────────────────────────────────────────────
"""
    prepare_df(raw) -> DataFrame

Keep the columns needed for ANOVA / SHAP and coerce every predictor to
`CategoricalArray`. `max_depth` is cast to String first so it is treated as
categorical rather than continuous.
"""
function prepare_df(raw::DataFrame)
    df = DataFrame(
        total_ec   = Float64.(raw.total_ec),
        name       = string.(raw.name),
        category   = string.(raw.category),
        row_idx    = string.(raw.row),
        col_idx    = string.(raw.col),
        slab_type  = string.(raw.slab_type),
        slab_sizer = string.(raw.slab_sizer),
        beam_sizer = string.(raw.beam_sizer),
        collinear  = string.(raw.collinear),
        max_depth  = string.(raw.max_depth),
    )
    for col in (:name, :category, :row_idx, :col_idx, :slab_type, :slab_sizer,
                :beam_sizer, :collinear, :max_depth)
        df[!, col] = categorical(df[!, col])
    end
    return df
end

"""
    TermSpec

Book-keeping record for one effect under test: a display `label`, a `kind`
tag used for plot colouring (`"layout"`, `"main"`, `"interaction"`), the
decision `Symbol`s involved (empty for the layout nuisance, one for a main
effect, two for a pairwise interaction), and the StatsModels `term`
object itself.
"""
struct TermSpec
    label::String
    kind::String
    syms::Vector{Symbol}
    term::Any
end

"""
    build_layout_nuisance_specs_flat(layout) -> Vector{TermSpec}

Original nuisance structure: one factor level per slab identity.
"""
function build_layout_nuisance_specs_flat(layout::Symbol)
    return TermSpec[
        TermSpec("layout ($(String(layout)))", "layout", Symbol[], Term(layout)),
    ]
end

"""
    build_layout_nuisance_specs_nested() -> Vector{TermSpec}

Hierarchical layout decomposition:
  - category
  - category × row_idx
  - category × col_idx
  - category × row_idx × col_idx
"""
function build_layout_nuisance_specs_nested()
    cat = Term(:category)
    row = Term(:row_idx)
    col = Term(:col_idx)
    return TermSpec[
        TermSpec("layout: category", "layout", Symbol[], cat),
        TermSpec("layout: category × row", "layout", Symbol[], cat & row),
        TermSpec("layout: category × col", "layout", Symbol[], cat & col),
        TermSpec("layout: category × row × col", "layout", Symbol[], cat & row & col),
    ]
end

"""
    build_term_specs(decisions; layout_mode, layout) -> Vector{TermSpec}

Layout nuisance block first, then decision main effects, then decision
pairwise interactions. This order keeps the decomposition table readable.
"""
function build_term_specs(decisions::Vector{Symbol};
                          layout_mode::Symbol = :flat,
                          layout::Symbol = LAYOUT)
    layout_specs = if layout_mode == :flat
        build_layout_nuisance_specs_flat(layout)
    elseif layout_mode == :nested
        build_layout_nuisance_specs_nested()
    else
        error("Unknown layout_mode=$(layout_mode). Use :flat or :nested.")
    end
    specs = copy(layout_specs)
    for d in decisions
        push!(specs, TermSpec(String(d), "main", [d], Term(d)))
    end
    for i in 1:length(decisions)-1, j in i+1:length(decisions)
        a, b = decisions[i], decisions[j]
        push!(specs, TermSpec("$(String(a)) × $(String(b))",
            "interaction", [a, b], Term(a) & Term(b)))
    end
    return specs
end

"""
    type_ii_specs(target, all_specs) -> (ref_specs, drop_specs)

Type II SS reduction rule, returned as the *reference* model (against which
`target` is tested) and the *drop* model (reference minus `target`).

- For a **main effect** A, the reference model first removes every pairwise
  interaction containing A (marginality), then the drop model removes A
  itself. SS(A) then isolates the main effect, with df equal to the main
  effect's own degrees of freedom.
- For an **interaction** or the **layout** nuisance, the reference is the
  full model and the drop model simply omits the target term.
"""
function type_ii_specs(target::TermSpec, all_specs::Vector{TermSpec})
    if target.kind in ("main", "layout_main")
        sym = target.syms[1]
        ref  = [s for s in all_specs
                if !(s.kind in ("interaction", "layout_interaction") && sym in s.syms)]
        drop = [s for s in ref if s.label != target.label]
    else
        ref  = all_specs
        drop = [s for s in all_specs if s.label != target.label]
    end
    return ref, drop
end

"""
    eta_squared_table(df; response, decisions, layout) -> (DataFrame, NamedTuple)

Fit the full fixed-effects model and, for every decision main effect, every
pairwise interaction and the layout nuisance factor, compute Type II SS,
degrees of freedom, F, p, η² and ω². Main effects are tested against a
reduced model that drops all higher-order interactions involving them, so
they are not hidden by their own interactions (the marginality rule).

Returns the term table sorted by η² descending plus a small named tuple
with overall fit diagnostics (n, R², MSE, SS_total, residual DoF).
"""
function eta_squared_table(df::DataFrame;
                           response::Symbol = :total_ec,
                           decisions::Vector{Symbol} = DECISIONS,
                           layout_mode::Symbol = :flat,
                           layout::Symbol = LAYOUT)

    all_specs = build_term_specs(decisions; layout_mode = layout_mode, layout = layout)
    full_rhs  = reduce(+, [s.term for s in all_specs])
    full_fm   = lm(FormulaTerm(Term(response), full_rhs), df)

    y         = df[!, response]
    ss_total  = sum((y .- mean(y)).^2)
    ssr_full  = deviance(full_fm)
    dof_full  = Int(dof_residual(full_fm))
    mse_full  = ssr_full / dof_full
    r2_full   = 1 - ssr_full / ss_total

    rows = DataFrame(term = String[], kind = String[],
                     df = Int[], SS = Float64[], MS = Float64[],
                     F = Float64[], p = Float64[],
                     eta2 = Float64[], omega2 = Float64[])

    # Cache fits keyed by a canonical label string so each intermediate model
    # is fit at most once across the whole loop.
    fit_cache = Dict{String, Any}()
    fit_cache["FULL"] = full_fm
    function fit_model(specs)
        key = join(sort([s.label for s in specs]), "|")
        get!(fit_cache, key) do
            lm(FormulaTerm(Term(response), reduce(+, [s.term for s in specs])), df)
        end
    end

    for spec in all_specs
        ref_specs, drop_specs = type_ii_specs(spec, all_specs)
        isempty(drop_specs) && continue

        ref_fm  = spec.kind == "main" ? fit_model(ref_specs) : full_fm
        drop_fm = fit_model(drop_specs)
        ssr_ref = deviance(ref_fm)
        dof_ref = Int(dof_residual(ref_fm))
        ssr_drp = deviance(drop_fm)
        dof_drp = Int(dof_residual(drop_fm))

        ss_term = ssr_drp - ssr_ref
        df_term = dof_drp - dof_ref
        if df_term <= 0
            @warn "Skipping $(spec.label) — zero df after reduction"
            continue
        end

        ms_term = ss_term / df_term
        # Use the full model's MSE as the error term so all F / p tests share a
        # common denominator even when the reference model differs between
        # main effects and interactions (standard Type II convention).
        F      = ms_term / mse_full
        p      = ccdf(FDist(df_term, dof_full), F)
        eta2   = ss_term / ss_total
        omega2 = (ss_term - df_term * mse_full) / (ss_total + mse_full)

        push!(rows, (term = spec.label, kind = spec.kind,
                     df = df_term, SS = ss_term, MS = ms_term,
                     F = F, p = p, eta2 = eta2, omega2 = omega2))
    end

    sort!(rows, :eta2, rev = true)
    return rows, (ssr_full = ssr_full, ss_total = ss_total,
                  mse_full = mse_full, r2 = r2_full,
                  dof_full = dof_full, n = nrow(df))
end

"""
    build_term_specs_per_category(decisions) -> Vector{TermSpec}

Term structure for within-category ANOVA:
  row_idx + col_idx + row_idx×col_idx + decision mains + decision pairwise.
"""
function build_term_specs_per_category(decisions::Vector{Symbol})
    specs = TermSpec[
        TermSpec("row", "layout_main", [:row_idx], Term(:row_idx)),
        TermSpec("col", "layout_main", [:col_idx], Term(:col_idx)),
        TermSpec("row × col", "layout_interaction", [:row_idx, :col_idx], Term(:row_idx) & Term(:col_idx)),
    ]
    for d in decisions
        push!(specs, TermSpec(String(d), "main", [d], Term(d)))
    end
    for i in 1:length(decisions)-1, j in i+1:length(decisions)
        a, b = decisions[i], decisions[j]
        push!(specs, TermSpec("$(String(a)) × $(String(b))",
            "interaction", [a, b], Term(a) & Term(b)))
    end
    return specs
end

"""
    eta_squared_table_per_category(df, category_value; response, decisions)

Run ANOVA on one category subset using:
  total_ec ~ row + col + row×col + decisions + decision pairwise interactions
"""
function eta_squared_table_per_category(df::DataFrame, category_value::AbstractString;
                                        response::Symbol = :total_ec,
                                        decisions::Vector{Symbol} = DECISIONS)
    sub = filter(row -> String(row.category) == category_value, df)
    nrow(sub) == 0 && error("No rows found for category='$category_value'")

    all_specs = build_term_specs_per_category(decisions)
    full_rhs  = reduce(+, [s.term for s in all_specs])
    full_fm   = lm(FormulaTerm(Term(response), full_rhs), sub)

    y         = sub[!, response]
    ss_total  = sum((y .- mean(y)).^2)
    ssr_full  = deviance(full_fm)
    dof_full  = Int(dof_residual(full_fm))
    mse_full  = ssr_full / dof_full
    r2_full   = 1 - ssr_full / ss_total

    rows = DataFrame(term = String[], kind = String[],
                     df = Int[], SS = Float64[], MS = Float64[],
                     F = Float64[], p = Float64[],
                     eta2 = Float64[], omega2 = Float64[])

    fit_cache = Dict{String, Any}()
    fit_cache["FULL"] = full_fm
    function fit_model(specs)
        key = join(sort([s.label for s in specs]), "|")
        get!(fit_cache, key) do
            lm(FormulaTerm(Term(response), reduce(+, [s.term for s in specs])), sub)
        end
    end

    for spec in all_specs
        ref_specs, drop_specs = type_ii_specs(spec, all_specs)
        isempty(drop_specs) && continue

        ref_fm  = spec.kind == "main" ? fit_model(ref_specs) : full_fm
        drop_fm = fit_model(drop_specs)
        ssr_ref = deviance(ref_fm)
        dof_ref = Int(dof_residual(ref_fm))
        ssr_drp = deviance(drop_fm)
        dof_drp = Int(dof_residual(drop_fm))

        ss_term = ssr_drp - ssr_ref
        df_term = dof_drp - dof_ref
        if df_term <= 0
            @warn "Skipping $(spec.label) in category=$category_value — zero df after reduction"
            continue
        end

        ms_term = ss_term / df_term
        F       = ms_term / mse_full
        p       = ccdf(FDist(df_term, dof_full), F)
        eta2    = ss_term / ss_total
        omega2  = (ss_term - df_term * mse_full) / (ss_total + mse_full)

        push!(rows, (term = spec.label, kind = spec.kind,
                     df = df_term, SS = ss_term, MS = ms_term,
                     F = F, p = p, eta2 = eta2, omega2 = omega2))
    end

    sort!(rows, :eta2, rev = true)
    return rows, (ssr_full = ssr_full, ss_total = ss_total,
                  mse_full = mse_full, r2 = r2_full,
                  dof_full = dof_full, n = nrow(sub))
end

"""
    nested_layout_breakdown(df; response, decisions) -> (DataFrame, NamedTuple)

Sequential (Type I) variance decomposition of layout hierarchy on top of the
decision model:

  baseline = decisions + decision pairwise interactions
  + category
  + category × row_idx
  + category × col_idx
  + category × row_idx × col_idx

This explicitly allocates the layout nuisance contribution in the requested
nested order, so the components are additive and interpretable.
"""
function nested_layout_breakdown(df::DataFrame;
                                 response::Symbol = :total_ec,
                                 decisions::Vector{Symbol} = DECISIONS)
    # Build decision-only term block (same decision structure as main ANOVA).
    decision_specs = TermSpec[]
    for d in decisions
        push!(decision_specs, TermSpec(String(d), "main", [d], Term(d)))
    end
    for i in 1:length(decisions)-1, j in i+1:length(decisions)
        a, b = decisions[i], decisions[j]
        push!(decision_specs, TermSpec("$(String(a)) × $(String(b))",
            "interaction", [a, b], Term(a) & Term(b)))
    end

    cat = Term(:category)
    row = Term(:row_idx)
    col = Term(:col_idx)
    nested_blocks = [
        ("layout: category", cat),
        ("layout: category × row", cat & row),
        ("layout: category × col", cat & col),
        ("layout: category × row × col", cat & row & col),
    ]

    y        = df[!, response]
    ss_total = sum((y .- mean(y)).^2)

    # Baseline (decisions only)
    base_rhs = reduce(+, [s.term for s in decision_specs])
    base_fm  = lm(FormulaTerm(Term(response), base_rhs), df)
    ssr_prev = deviance(base_fm)
    dof_prev = Int(dof_residual(base_fm))

    # Full nested model (for common error denominator).
    all_terms = vcat([s.term for s in decision_specs], [b[2] for b in nested_blocks])
    full_fm   = lm(FormulaTerm(Term(response), reduce(+, all_terms)), df)
    ssr_full  = deviance(full_fm)
    dof_full  = Int(dof_residual(full_fm))
    mse_full  = ssr_full / dof_full
    r2_full   = 1 - ssr_full / ss_total

    rows = DataFrame(term = String[], kind = String[],
                     df = Int[], SS = Float64[], MS = Float64[],
                     F = Float64[], p = Float64[],
                     eta2 = Float64[], omega2 = Float64[])

    current_terms = [s.term for s in decision_specs]
    for (label, term) in nested_blocks
        push!(current_terms, term)
        cur_fm  = lm(FormulaTerm(Term(response), reduce(+, current_terms)), df)
        ssr_cur = deviance(cur_fm)
        dof_cur = Int(dof_residual(cur_fm))

        ss_term = ssr_prev - ssr_cur
        df_term = dof_prev - dof_cur
        ms_term = ss_term / df_term
        F       = ms_term / mse_full
        p       = ccdf(FDist(df_term, dof_full), F)
        eta2    = ss_term / ss_total
        omega2  = (ss_term - df_term * mse_full) / (ss_total + mse_full)

        push!(rows, (term = label, kind = "layout",
                     df = df_term, SS = ss_term, MS = ms_term,
                     F = F, p = p, eta2 = eta2, omega2 = omega2))

        ssr_prev = ssr_cur
        dof_prev = dof_cur
    end

    # Convenience total for layout hierarchy sum.
    layout_total_ss = sum(rows.SS)
    layout_total_df = sum(rows.df)
    layout_total_ms = layout_total_ss / layout_total_df
    layout_total_F  = layout_total_ms / mse_full
    layout_total_p  = ccdf(FDist(layout_total_df, dof_full), layout_total_F)
    layout_total_eta2 = layout_total_ss / ss_total
    layout_total_omega2 = (layout_total_ss - layout_total_df * mse_full) / (ss_total + mse_full)
    push!(rows, (term = "layout: total nested", kind = "layout",
                 df = layout_total_df, SS = layout_total_ss, MS = layout_total_ms,
                 F = layout_total_F, p = layout_total_p,
                 eta2 = layout_total_eta2, omega2 = layout_total_omega2))

    return rows, (ss_total = ss_total, ssr_full = ssr_full, mse_full = mse_full,
                  dof_full = dof_full, r2 = r2_full, n = nrow(df))
end

"""
    plot_eta2_bar(rows, fit, label, outpath)

Horizontal bar chart of η² per term. Layout is drawn as a muted reference
bar at the top; decision main effects and pairwise interactions are colour-
coded so reviewers can see at a glance whether variance lives in linear
factors or in interactions. Bars are labelled with the percentage.
"""
function plot_eta2_bar(rows::DataFrame, fit, label::String, outpath::String)
    n_rows = nrow(rows)
    colour = Dict(
        "layout"      => 色[:charcoalgrey],
        "layout_main" => 色[:charcoalgrey],
        "layout_interaction" => colorant"#6b6b6b",
        "main"        => 色[:ceruleanblue],
        "interaction" => 色[:magenta],
    )
    bar_colours = [colour[k] for k in rows.kind]
    bar_vals    = rows.eta2 .* 100

    fig = Figure(size = (780, 60 + 28 * n_rows), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "η² (% of total EC variance)",
        title  = "Variance decomposition — $label   (n = $(fit.n),  R² = $(round(fit.r2, digits=3)))",
        yticks = (1:n_rows, reverse(rows.term)),
        yticklabelsize = SMALLFONTSIZE,
    )

    barplot!(ax, n_rows:-1:1, bar_vals,
        direction = :x,
        color = bar_colours,
        strokewidth = 0,
    )

    # Percentage labels just to the right of each bar end.
    for (i, v) in enumerate(bar_vals)
        text!(ax, "$(round(v, digits=2))%",
            position = (v, n_rows - i + 1),
            align = (:left, :center),
            offset = (4, 0),
            fontsize = SMALLFONTSIZE)
    end

    # Manual legend (three swatches).
    elems = [PolyElement(color = colour[k]) for k in ("layout", "main", "interaction")]
    Legend(fig[1, 2], elems,
        ["layout (108 slabs)", "decision main effect", "pairwise interaction"],
        framevisible = false, labelsize = SMALLFONTSIZE)

    xlims!(ax, (0, max(maximum(bar_vals) * 1.15, 5.0)))
    save(outpath, fig)
    return fig
end

"""
    plot_eta2_comparison(rows_a, rows_b, labels, outpath)

Grouped horizontal bar chart comparing η² for each shared term across the
two reviewer conditions. Terms are drawn in the union of both tables, sorted
by the larger of the two η² values so the dominant factors anchor the top.
"""
function plot_eta2_comparison(rows_a::DataFrame, rows_b::DataFrame,
                              labels::NTuple{2,String}, outpath::String)
    terms = unique(vcat(rows_a.term, rows_b.term))
    lookup(r, t) = (idx = findfirst(==(t), r.term);
                    idx === nothing ? 0.0 : r.eta2[idx] * 100)

    vals_a = [lookup(rows_a, t) for t in terms]
    vals_b = [lookup(rows_b, t) for t in terms]
    order  = sortperm(max.(vals_a, vals_b), rev = true)
    terms  = terms[order]
    vals_a = vals_a[order]
    vals_b = vals_b[order]
    n_rows = length(terms)

    fig = Figure(size = (820, 80 + 30 * n_rows), fontsize = FONTSIZE)
    ax  = Axis(fig[1, 1],
        xlabel = "η² (% of total EC variance)",
        title  = "Variance decomposition — condition comparison",
        yticks = (1:n_rows, reverse(terms)),
        yticklabelsize = SMALLFONTSIZE,
    )

    bar_w = 0.36
    ys    = n_rows:-1:1
    barplot!(ax, ys .+ bar_w/2, vals_a, direction = :x,
        color = 色[:ceruleanblue], strokewidth = 0, width = bar_w,
        label = labels[1])
    barplot!(ax, ys .- bar_w/2, vals_b, direction = :x,
        color = 色[:magenta], strokewidth = 0, width = bar_w,
        label = labels[2])

    axislegend(ax, position = :rb, labelsize = SMALLFONTSIZE,
               framevisible = false)
    xlims!(ax, (0, max(maximum(vals_a), maximum(vals_b)) * 1.15 + 1))
    save(outpath, fig)
    return fig
end

# ── Python companion: written once under SAVE_ROOT for both conditions ───────
const SHAP_SCRIPT_CONTENT = raw"""
#!/usr/bin/env python3
'''SHAP analysis on the master design matrix exported by
`sensitivity_variance_decomposition.jl`.

Run separately for each condition from that condition's subfolder, e.g.

    cd SlabDesignFactors/results/sensitivity_variance_decomposition/yesdeflection_yesslabmin
    python ../shap_analysis.py master_design_matrix.csv

Requires `lightgbm`, `shap`, `pandas`, `matplotlib`. Outputs:
  - shap_summary_bar.pdf       global importance (mean |SHAP|)
  - shap_summary_beeswarm.pdf  per-sample SHAP with direction of effect
  - shap_values.csv            raw SHAP value matrix (rows × features)
'''
import sys, os
import pandas as pd
import numpy as np
import lightgbm as lgb
import shap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_path = sys.argv[1] if len(sys.argv) > 1 else "master_design_matrix.csv"
outdir   = os.path.dirname(os.path.abspath(csv_path)) or "."

df = pd.read_csv(csv_path)
y  = df["total_ec"].values
X  = df.drop(columns=["total_ec"])

# one-hot encode categorical columns; LightGBM can also handle them natively
# via `categorical_feature`, but one-hot keeps SHAP interpretation simple.
X_enc = pd.get_dummies(X, drop_first=False)

model = lgb.LGBMRegressor(
    n_estimators=800,
    learning_rate=0.05,
    num_leaves=63,
    min_data_in_leaf=8,
    feature_fraction=0.9,
    bagging_fraction=0.9,
    bagging_freq=3,
    random_state=42,
)
model.fit(X_enc, y)

print(f"Train R²: {model.score(X_enc, y):.3f}")

explainer = shap.TreeExplainer(model)
shap_vals = explainer.shap_values(X_enc)

# Save raw SHAP matrix
shap_df = pd.DataFrame(shap_vals, columns=X_enc.columns)
shap_df.to_csv(os.path.join(outdir, "shap_values.csv"), index=False)

# Bar plot — global importance
plt.figure()
shap.summary_plot(shap_vals, X_enc, plot_type="bar", show=False, max_display=25)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "shap_summary_bar.pdf"))
plt.close()

# Beeswarm — per-sample SHAP distribution
plt.figure()
shap.summary_plot(shap_vals, X_enc, show=False, max_display=25)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "shap_summary_beeswarm.pdf"))
plt.close()

print("✓ SHAP outputs written to", outdir)
"""

function write_shap_companion(path::String)
    open(path, "w") do io
        write(io, SHAP_SCRIPT_CONTENT)
    end
end

# ── main loop ─────────────────────────────────────────────────────────────────
condition_results_flat   = Dict{String, DataFrame}()
condition_results_nested = Dict{String, DataFrame}()
condition_results_nested_breakdown = Dict{String, DataFrame}()
condition_results_per_category = Dict{String, Dict{String, DataFrame}}()
condition_labels         = Dict{String, String}()

for cond in CONDITIONS
    println("\n════════════════════════════════════════════")
    println("  Condition: $(cond.label)")
    println("  Source   : $(cond.base)")
    println("════════════════════════════════════════════")

    out_dir = joinpath(SAVE_ROOT, cond.tag)
    mkpath(out_dir)

    raw = assemble_data(cond.base)
    raw = filter(row -> row.area > 0 && !isnan(row.steel_norm), raw)
    df  = prepare_df(raw)
    println("Prepared $(nrow(df)) rows  ·  $(length(levels(df.name))) distinct layouts")

    # ANOVA A: original flat layout nuisance (preserves prior result)
    rows_flat, fit_flat = eta_squared_table(df; layout_mode = :flat, layout = LAYOUT)
    condition_results_flat[cond.tag] = rows_flat
    condition_labels[cond.tag]  = cond.label

    eta_csv_flat = joinpath(out_dir, "anova_eta_squared.csv")
    CSV.write(eta_csv_flat, rows_flat)
    println("✓ η² table (flat layout) → $eta_csv_flat")

    eta_pdf_flat = joinpath(out_dir, "anova_eta_squared.pdf")
    plot_eta2_bar(rows_flat, fit_flat, cond.label, eta_pdf_flat)
    println("✓ η² bar plot (flat layout) → $eta_pdf_flat")

    # ANOVA B: supplemental nested layout nuisance decomposition
    rows_nested, fit_nested = eta_squared_table(df; layout_mode = :nested)
    condition_results_nested[cond.tag] = rows_nested

    eta_csv_nested = joinpath(out_dir, "anova_eta_squared_nested_layout.csv")
    CSV.write(eta_csv_nested, rows_nested)
    println("✓ η² table (nested layout) → $eta_csv_nested")

    eta_pdf_nested = joinpath(out_dir, "anova_eta_squared_nested_layout.pdf")
    plot_eta2_bar(rows_nested, fit_nested, cond.label * "  [nested layout]", eta_pdf_nested)
    println("✓ η² bar plot (nested layout) → $eta_pdf_nested")

    # Supplemental sequential nested-layout breakdown (explicit additive split).
    rows_nested_break, fit_nested_break = nested_layout_breakdown(df)
    condition_results_nested_breakdown[cond.tag] = rows_nested_break
    nested_break_csv = joinpath(out_dir, "anova_layout_nested_breakdown.csv")
    CSV.write(nested_break_csv, rows_nested_break)
    println("✓ Layout nested breakdown (sequential) → $nested_break_csv")

    nested_break_pdf = joinpath(out_dir, "anova_layout_nested_breakdown.pdf")
    plot_eta2_bar(rows_nested_break, fit_nested_break,
        cond.label * "  [layout nested breakdown]", nested_break_pdf)
    println("✓ Layout nested breakdown plot → $nested_break_pdf")

    # ANOVA C: per-category ANOVAs (grid / nova / topology), each with
    # row + col + row×col + decisions + decision pairwise interactions.
    categories = sort(collect(String.(levels(df.category))))
    per_cat = Dict{String, DataFrame}()
    per_cat_long = DataFrame(
        category = String[], n = Int[], r2 = Float64[],
        term = String[], kind = String[],
        df = Int[], SS = Float64[], MS = Float64[], F = Float64[],
        p = Float64[], eta2 = Float64[], omega2 = Float64[],
    )

    for cat in categories
        rows_cat, fit_cat = eta_squared_table_per_category(df, cat)
        per_cat[cat] = rows_cat

        cat_csv = joinpath(out_dir, "anova_eta_squared_category_$(cat).csv")
        CSV.write(cat_csv, rows_cat)
        println("✓ η² table (category=$(cat)) → $cat_csv")

        cat_pdf = joinpath(out_dir, "anova_eta_squared_category_$(cat).pdf")
        plot_eta2_bar(rows_cat, fit_cat, cond.label * "  [category=$(cat)]", cat_pdf)
        println("✓ η² bar plot (category=$(cat)) → $cat_pdf")

        cat_long = copy(rows_cat)
        cat_long[!, :category] .= cat
        cat_long[!, :n] .= fit_cat.n
        cat_long[!, :r2] .= fit_cat.r2
        append!(per_cat_long, cat_long[:, [:category, :n, :r2, :term, :kind,
                                           :df, :SS, :MS, :F, :p, :eta2, :omega2]])
    end
    condition_results_per_category[cond.tag] = per_cat

    long_csv = joinpath(out_dir, "anova_eta_squared_by_category_long.csv")
    CSV.write(long_csv, per_cat_long)
    println("✓ Per-category long table → $long_csv")

    terms_union = sort(unique(per_cat_long.term))
    wide = DataFrame(term = terms_union)
    for cat in categories
        tbl = per_cat[cat]
        function _lookup_term(term, col)
            idx = findfirst(==(term), tbl.term)
            return idx === nothing ? NaN : tbl[idx, col]
        end
        wide[!, Symbol("eta2_" * cat)] = [_lookup_term(t, :eta2) for t in terms_union]
        wide[!, Symbol("omega2_" * cat)] = [_lookup_term(t, :omega2) for t in terms_union]
    end
    wide_csv = joinpath(out_dir, "anova_eta_squared_by_category_wide.csv")
    CSV.write(wide_csv, wide)
    println("✓ Per-category wide table → $wide_csv")

    # Fit diagnostics
    open(joinpath(out_dir, "fit_summary.txt"), "w") do io
        println(io, "Condition       : $(cond.label)")
        println(io, "Source directory: $(cond.base)")
        println(io, "n (observations): $(fit_flat.n)")
        println(io, "SS_total        : $(fit_flat.ss_total)")
        println(io, "")
        println(io, "Model A (flat layout nuisance)")
        println(io, "  layout term   : $(LAYOUT)")
        println(io, "  SS_residual   : $(fit_flat.ssr_full)")
        println(io, "  MSE           : $(fit_flat.mse_full)")
        println(io, "  R²            : $(fit_flat.r2)")
        println(io, "  DoF residual  : $(fit_flat.dof_full)")
        println(io, "")
        println(io, "Model B (nested layout nuisance)")
        println(io, "  layout terms  : category + category×row + category×col + category×row×col")
        println(io, "  SS_residual   : $(fit_nested.ssr_full)")
        println(io, "  MSE           : $(fit_nested.mse_full)")
        println(io, "  R²            : $(fit_nested.r2)")
        println(io, "  DoF residual  : $(fit_nested.dof_full)")
        println(io, "Response        : total_ec (kgCO₂e/m²)")
        println(io, "Primary result  : anova_eta_squared.csv (flat layout)")
        println(io, "Supplement      : anova_eta_squared_nested_layout.csv")
        println(io, "Nested breakdown: anova_layout_nested_breakdown.csv")
        println(io, "Per-category    : anova_eta_squared_by_category_long.csv")
        println(io, "Per-category    : anova_eta_squared_by_category_wide.csv")
        decisions_str = join(DECISIONS, ", ")
        println(io, "Decision factors: $decisions_str")
    end

    # Export master design matrix for downstream SHAP.
    master = DataFrame(
        total_ec   = df.total_ec,
        category   = raw.category,
        layout     = df.name,
        row        = raw.row,
        col        = raw.col,
        slab_type  = df.slab_type,
        slab_sizer = df.slab_sizer,
        beam_sizer = df.beam_sizer,
        collinear  = df.collinear,
        max_depth  = df.max_depth,
    )
    master_csv = joinpath(out_dir, "master_design_matrix.csv")
    CSV.write(master_csv, master)
    println("✓ Master design matrix → $master_csv  ($(nrow(master)) rows)")
end

# Cross-condition comparison
tag_a, tag_b = CONDITIONS[1].tag, CONDITIONS[2].tag
rows_a_flat, rows_b_flat = condition_results_flat[tag_a], condition_results_flat[tag_b]
rows_a_nested, rows_b_nested = condition_results_nested[tag_a], condition_results_nested[tag_b]
rows_a_nested_break = condition_results_nested_breakdown[tag_a]
rows_b_nested_break = condition_results_nested_breakdown[tag_b]

combined_flat = outerjoin(
    rename(rows_a_flat[:, [:term, :kind, :eta2, :omega2]],
           :eta2 => :eta2_yesslabmin, :omega2 => :omega2_yesslabmin),
    rename(rows_b_flat[:, [:term, :eta2, :omega2]],
           :eta2 => :eta2_noslabmin,  :omega2 => :omega2_noslabmin),
    on = :term,
)
sort!(combined_flat,
      [order(:eta2_yesslabmin, rev = true, by = x -> ismissing(x) ? 0.0 : x)])
CSV.write(joinpath(SAVE_ROOT, "eta_squared_comparison.csv"), combined_flat)

plot_eta2_comparison(rows_a_flat, rows_b_flat,
    (condition_labels[tag_a], condition_labels[tag_b]),
    joinpath(SAVE_ROOT, "eta_squared_comparison.pdf"))

combined_nested = outerjoin(
    rename(rows_a_nested[:, [:term, :kind, :eta2, :omega2]],
           :eta2 => :eta2_yesslabmin, :omega2 => :omega2_yesslabmin),
    rename(rows_b_nested[:, [:term, :eta2, :omega2]],
           :eta2 => :eta2_noslabmin,  :omega2 => :omega2_noslabmin),
    on = :term,
)
sort!(combined_nested,
      [order(:eta2_yesslabmin, rev = true, by = x -> ismissing(x) ? 0.0 : x)])
CSV.write(joinpath(SAVE_ROOT, "eta_squared_comparison_nested_layout.csv"), combined_nested)

plot_eta2_comparison(rows_a_nested, rows_b_nested,
    (condition_labels[tag_a], condition_labels[tag_b]),
    joinpath(SAVE_ROOT, "eta_squared_comparison_nested_layout.pdf"))

combined_nested_break = outerjoin(
    rename(rows_a_nested_break[:, [:term, :kind, :eta2, :omega2]],
           :eta2 => :eta2_yesslabmin, :omega2 => :omega2_yesslabmin),
    rename(rows_b_nested_break[:, [:term, :eta2, :omega2]],
           :eta2 => :eta2_noslabmin,  :omega2 => :omega2_noslabmin),
    on = :term,
)
sort!(combined_nested_break,
      [order(:eta2_yesslabmin, rev = true, by = x -> ismissing(x) ? 0.0 : x)])
CSV.write(joinpath(SAVE_ROOT, "eta_squared_comparison_nested_breakdown.csv"), combined_nested_break)

plot_eta2_comparison(rows_a_nested_break, rows_b_nested_break,
    (condition_labels[tag_a], condition_labels[tag_b]),
    joinpath(SAVE_ROOT, "eta_squared_comparison_nested_breakdown.pdf"))

# Per-category, cross-condition wide table (6 η² columns = 3 categories × 2 conditions).
all_terms_cat = String[]
for cond in CONDITIONS
    per_cat = condition_results_per_category[cond.tag]
    for tbl in values(per_cat)
        append!(all_terms_cat, tbl.term)
    end
end
all_terms_cat = sort(unique(all_terms_cat))
wide_cat = DataFrame(term = all_terms_cat)

for cond in CONDITIONS
    per_cat = condition_results_per_category[cond.tag]
    for cat in sort(collect(keys(per_cat)))
        tbl = per_cat[cat]
        function _lookup_term(term, col)
            idx = findfirst(==(term), tbl.term)
            return idx === nothing ? NaN : tbl[idx, col]
        end
        wide_cat[!, Symbol("eta2_" * cond.tag * "_" * cat)] =
            [_lookup_term(t, :eta2) for t in all_terms_cat]
        wide_cat[!, Symbol("omega2_" * cond.tag * "_" * cat)] =
            [_lookup_term(t, :omega2) for t in all_terms_cat]
    end
end
CSV.write(joinpath(SAVE_ROOT, "eta_squared_by_category_across_conditions.csv"), wide_cat)

# Drop the Python SHAP companion once, at the top level.
write_shap_companion(joinpath(SAVE_ROOT, "shap_analysis.py"))

println("\n════════════════════════════════════════════")
println("  Done. Outputs under:")
println("    $(SAVE_ROOT)")
println("  For SHAP:  cd <condition>/ && python ../shap_analysis.py master_design_matrix.csv")
println("════════════════════════════════════════════")
