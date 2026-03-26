"""
Verify that consolidate_loads preserves global equilibrium.

For each slab type (isotropic, uniaxial, orth_biaxial), run analyze_slab on
a simple r1c1 geometry and check vertical equilibrium: sum(applied) + sum(reactions) ≈ 0.

`sum(slab_params.load_volumes)` is **not** equal to total slab volume (area × depth) on
multi-cycle topologies: each cycle appends strip volumes, so interior contributions
appear twice in that global sum. After `consolidate_loads`, duplicate beam stations
are merged with the `/2` rule, so **post/pre ≈ 0.5** is expected here — not a bug.

The authoritative total for the solved model is the **post-consolidation** sum of
`model.loads`, which matches reactions (equilibrium).
"""

using Test

main_path = joinpath(@__DIR__, "..", "..", "Geometries", "topology")

@testset "consolidate_loads equilibrium" begin

    for slab_type in [:isotropic, :uniaxial, :orth_biaxial]
        @testset "$slab_type" begin
            json_path = joinpath(main_path, "r1c1.json")
            geometry_dict = geometry_dict_from_json_path(json_path)
            geom, _ = Base.invokelatest(generate_from_json, geometry_dict; plot=false, drawn=false)

            vector_1d = slab_type == :isotropic ? [0.0, 0.0] : [1.0, 0.0]

            slab_params = SlabAnalysisParams(
                geom,
                slab_name     = "equilibrium_test",
                slab_type     = slab_type,
                vector_1d     = vector_1d,
                slab_sizer    = :uniform,
                spacing       = 0.1,
                plot_analysis = false,
                fix_param     = true,
                slab_units    = :m,
            )

            slab_params = analyze_slab(slab_params)

            model = slab_params.model
            n_loads = length(model.loads)
            @test n_loads > 0

            # Sum of applied vertical point-load values (should be negative = downward)
            applied_vertical = sum(load.value[3] for load in model.loads)

            # Sum of vertical reactions at fixed (ground) nodes
            # In a solved model, reactions are stored at fixed-DOF nodes.
            # Total vertical reaction should equal -applied_vertical (equilibrium).
            # We check by summing reaction forces at ground nodes.
            reaction_vertical = 0.0
            for node in model.nodes
                if node.id == :fixed
                    # Reaction = external force at constrained DOFs after solve
                    reaction_vertical += node.reaction[3]
                end
            end

            println("\n[$slab_type] Applied vertical load: $(round(applied_vertical, digits=6))")
            println("[$slab_type] Reaction vertical:     $(round(reaction_vertical, digits=6))")
            println("[$slab_type] Sum (should be ~0):     $(round(applied_vertical + reaction_vertical, digits=6))")

            # Equilibrium: applied + reactions ≈ 0 (strict — must hold regardless of /2 bookkeeping)
            @test abs(applied_vertical + reaction_vertical) < 1e-6 * max(1.0, abs(applied_vertical))

            # Now check against expected load from geometry.
            # For density-based slab DL: each strip contributes volume = distance * spacing * slab_depth
            # Total load = sum(load_volumes) = sum(load_areas * slab_depth_per_cell)
            # But after consolidation the model loads are the post-processed values.
            # Instead, verify that slab_params.area > 0 and reactions are sensible.

            slab_area = slab_params.area  # m²
            @test slab_area > 0

            # The applied load magnitudes should be proportional to area.
            # For isotropic uniform slab, total applied ≈ -slab_area * mean_slab_depth
            # (since each load.value[3] = volume = area_strip * depth, and depth is the "w" passed in)
            mean_depth = isempty(slab_params.slab_depths) ? 0.0 : sum(slab_params.slab_depths) / length(slab_params.slab_depths)

            if slab_type == :isotropic && slab_params.slab_sizer == :uniform
                naive_volume_sum = -slab_area * mean_depth  # matches sum(load_volumes) when strips partition the slab once per cycle list
                println("[$slab_type] Slab area: $(round(slab_area, digits=2)) m²")
                println("[$slab_type] Mean depth: $(round(mean_depth, digits=4)) m")
                println("[$slab_type] Naive area×depth (upper bound): $(round(naive_volume_sum, digits=6))")
                println("[$slab_type] Actual consolidated load:     $(round(applied_vertical, digits=6))")
                println("[$slab_type] Ratio consolidated / (area×depth): $(round(applied_vertical / naive_volume_sum, digits=4))")
            end

            # Also directly verify: sum of load_areas * corresponding depths ≈ total applied
            # load_areas and load_volumes are the PRE-consolidation strip values
            if !isempty(slab_params.load_volumes)
                pre_consolidation_total = -sum(slab_params.load_volumes)
                post_consolidation_total = applied_vertical
                ratio_pre_post = post_consolidation_total / pre_consolidation_total
                println("[$slab_type] Pre-consolidation total:  $(round(pre_consolidation_total, digits=6))")
                println("[$slab_type] Post-consolidation total: $(round(post_consolidation_total, digits=6))")
                println("[$slab_type] Ratio post/pre:           $(round(ratio_pre_post, digits=4))")

                # On multi-cycle geometries, sum(load_volumes) double-counts strip volumes
                # across cycles; consolidate_loads merges stations and /2 averages paired
                # contributions — expect post/pre ≈ 0.5 here (r1c1), not 1.0.
                @test 0.45 < ratio_pre_post < 0.55
            end
        end
    end
end
