"""
    FULL_SWEEP_* constants

Single source of truth for **topology** studies that share one load/factor/strip
policy: `executable_experiments.run_max_depths`, `strip_resolution`,
`nlp_solver_comparison`, `material_scenario_mc`, and
`executable_analyze_sharded` (full sweep).

Validation / constrained-inventory studies use building-specific loads instead;
do not reuse these there.

Notes
-----
- `SlabSizingParams.deflection_limit` is left at its struct default (`true`) for
  these sweeps.
- `slab_dead_load` is left at `0.0` so slab self-weight follows material density
  when applicable.
- `composite_action`, `collinear`, and `deflection_reduction_factor` match the
  analyze → `optimal_beamsizer` → `postprocess_slab` pipeline exercised by
  `SlabDesignFactors/test/run.jl` (integration tests) and the cluster strip-resolution /
  NLP-comparison studies: composite stiffness, staged load decomposition in
  `get_scaled_model` when applicable, and L/360 + L/240 staged checks in the sizer
  and postprocessor.
"""
const FULL_SWEEP_LIVE_LOAD = psf_to_ksi(50)
const FULL_SWEEP_SUPERIMPOSED_DEAD_LOAD = psf_to_ksi(15)
const FULL_SWEEP_LIVE_FACTOR = 1.6
const FULL_SWEEP_DEAD_FACTOR = 1.2
const FULL_SWEEP_SERVICEABILITY_LIM = 360
const FULL_SWEEP_STRIP_SPACING = 0.01  # `SlabAnalysisParams.spacing` [m]

const FULL_SWEEP_COMPOSITE_ACTION = true
const FULL_SWEEP_COLLINEAR = true
const FULL_SWEEP_DEFLECTION_REDUCTION_FACTOR = 1.0
