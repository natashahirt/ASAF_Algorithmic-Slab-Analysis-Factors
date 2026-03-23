include("SlabDesignFactors/core/conversions.jl")
include("SlabDesignFactors/core/constants.jl")
include("SlabDesignFactors/core/materials.jl")
include("SlabDesignFactors/core/structs.jl")
include("TributaryAreas/slab_analysis/output.jl")
r = SlabOptimResults(
    slab_name = "t",
    δ_slab_dead = [0.01],
    δ_beam_dead = [0.02],
    δ_sdl = [0.03],
    δ_live = [0.04],
    δ_total = [0.05],
    Δ_limit_live = [1.0],
    Δ_limit_total = [1.0],
    δ_live_ok = [true],
    δ_total_ok = [false],
    n_L360_fail = 0,
    n_L240_fail = 1,
    i_L240_fail = [1],
    composite_action = true,
    staged_converged = false,
    staged_n_violations = 2,
    max_bay_span = 120.0,
    global_δ_ok = false,
    max_util_M = 0.9,
)
df = create_results_dataframe([r], false)
@assert "staged_converged" in names(df)
@assert df.staged_n_violations[1] == 2
@assert df.n_L240_fail[1] == 1
println("ok")
