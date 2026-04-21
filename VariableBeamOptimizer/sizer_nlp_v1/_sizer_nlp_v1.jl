"""
    sizer_nlp_v1

Validation-only port of the November 2024 bare-steel two-stage beam
sizing workflow as it existed at the original git state.

# Authoritative sources

  - `TributaryAreas @ a8dd2b98`  (2024-11-19) — the complete sizer
      **driver** as it ran on the cluster for the Nov 24 2024 sweep:
      `slab_analysis/size_beams.jl` (`iterate_discrete_continuous`,
      `process_discrete_beams`, `process_continuous_beams`) and
      `slab_analysis/size_beams_utils.jl` (`binary_search_sections`,
      `get_deflection_constraint` 3-arg form, `check_uniqueness!`).
      The monorepo's `SlabDesignFactors @ 2ede3f0` submodule pointer
      at Nov 24 2024 pins TributaryAreas exactly at this commit.
  - `VariableBeamOptimizer @ c3a506f`  (2024-11-13) — the NLP kernel:
      `optim/optimize_I_symm.jl` (`generate_broadcast`,
      `generate_objective` with soft `(δ-δ_max)^2 / (δ-δ_max)^4`
      penalty, and strength-only `inequality_constraints_I_symm`) and
      `optim/utils.jl` (`get_element_deflection`, `get_frameoptparams`).
      Also pinned by `SlabDesignFactors @ 2ede3f0`.

Together these two define, exactly, the Nov 2024 regime that generated
`SlabDesignFactors/results/remote_results_*deflection_*slabmin/topology.csv`.

# What this module contributes

The NLP kernel (`generate_broadcast`, `generate_objective`,
`inequality_constraints_I_symm`, `get_element_deflection`,
`get_frameoptparams`) already lives under
`VariableBeamOptimizer/rolled_steel/` byte-for-byte with c3a506f.
This module re-homes the TributaryAreas-era driver so it can be called
against the current `Asap`/`Nonconvex` versions:

  - `process_discrete_beams_nov2024` — ports the
    `binary_search_sections`-based discrete W-catalog sizer from
    a8dd2b98 (as opposed to the monorepo's current
    `sequential_search_sections` variant).
  - `process_continuous_beams_nov2024` — ports the continuous NLP
    driver from a8dd2b98, with:
      * `W4X13` as the floor for the NLP variable bounds when
        `minimum_continuous=true` (NOT `W6X8.5` — that was a later
        drift), and also as the cold-start section if no warm start
        is supplied.
      * `W43X335` as the ceiling.
      * `NLoptAlg(:LD_MMA)` with
        `xtol_rel = ftol_abs = ftol_rel = 1e-2, xtol_abs = 1e-8`.
      * Hard `get_deflection_constraint` pushed on top of the soft
        penalty in `generate_objective`.
  - `size_bare_nlp!` — a facade that reproduces a8dd2b98's
    `iterate_discrete_continuous`: run `process_discrete_beams_nov2024`
    first, harvest its minimizers, then warm-start
    `process_continuous_beams_nov2024` from them. Cold-starting
    without the discrete stage leaves MMA at a local basin ~35 %
    lighter than the warm-started global.

# Composite-aware-hook suppression

Three settings are pinned so the modern composite pipeline's hooks
stay out of the way:

  - `params.composite_action = false` and
    `params.deflection_reduction_factor = 1.0` (the facade overrides
    whatever the caller passed). This forces today's
    `get_deflection_constraint` down its bare-steel, non-composite
    fallback, which is numerically equivalent to a8dd2b98's simpler
    3-arg formula (`abs(minimum(δ_FE)) ≤ L/serviceability_lim`): the
    modern fallback passes `dead_load=0.0` to `get_element_deflection`
    and adds `δ_sw = 5·w_sw·L⁴/(384·E·Ix)` analytically, which for a
    prismatic simply-supported beam is algebraically identical to
    a8dd2b98's `dead_load=ρ` FE inclusion.
  - `params.serviceability_lim = 360` (overrideable).
  - No staged / iterative deflection loop: Nov 2024 did not have one.

# Not for production

Only the comparison script
`SlabDesignFactors/scripts/compare_nlp_vs_composite.jl` and its
reproduction sanity test call into this namespace.
"""

include("get_deflection_constraint_nov2024.jl")
include("process_discrete_beams_nov2024.jl")
include("process_continuous_beams_nov2024.jl")
include("size_bare_nlp.jl")
