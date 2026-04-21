#!/usr/bin/env julia
# Diagnostic: why does postprocess_slab report max_δ_total = 3.29" when Nov 2024
# sections were sized to bare-steel δ(total) ≤ L/360 ≈ 0.55"?
#
# Strategy: run size_bare_nlp! for BaU, then manually re-create the staged solve
# steps (A and B) on beam 1 and print every intermediate quantity so we can tell
# whether the load magnitudes, the Ix values, or something else is off.

include(joinpath(@__DIR__, "..", "scripts", "_scripts.jl"))

using JSON, Asap

const NAME = "r6c4"
const GEOM = joinpath(@__DIR__, "..", "..", "Geometries", "topology", "$NAME.json")

# Load geometry (same idiom as comparison script)
raw_json = JSON.parse(JSON.parse(replace(read(GEOM, String), "\\n" => ""); dicttype=Dict))
geom_dict = raw_json isa Dict ? raw_json : Dict(pairs(raw_json))
geom, _ = Base.invokelatest(generate_from_json, geom_dict; plot=false, drawn=false)

slab = SlabAnalysisParams(geom; slab_name=NAME, slab_type=:isotropic,
                          vector_1d=[0.0,0.0], slab_sizer=:uniform,
                          spacing=0.1, plot_analysis=false, fix_param=true, slab_units=:m)

sp  = SlabSizingParams(live_load=psf_to_ksi(50), superimposed_dead_load=psf_to_ksi(15),
                       slab_dead_load=0.0, live_factor=1.6, dead_factor=1.2,
                       beam_sizer=:continuous, nlp_solver=:MMA, max_depth=25.0,
                       beam_units=:in, serviceability_lim=360, collinear=false,
                       minimum_continuous=true, n_max_sections=0,
                       composite_action=false)

slab = analyze_slab(slab)
println("\n── Step 1: Nov 2024 NLP sizing ──")
slab, sp = size_bare_nlp!(slab, sp)

# Capture NLP solution and bare section properties for beam 1
nlp_mins = deepcopy(sp.minimizers)
nlp_ids  = deepcopy(sp.ids)
n_beams  = length(sp.model.elements[:beam])
L1_in    = sp.model.elements[:beam][1].length
sec1     = I_symm(nlp_mins[1]...)
bare_Ix1 = sec1.Ix

println("\n── Step 2: bare-steel δ(total unfactored) on beam 1 ──")
# Put sized sections onto the model's beam elements, apply unfactored (:all)
# loads, solve, measure δ. This is the closest analog to what the Nov 2024
# NLP constraint enforces (δ_bare ≤ L/360 with unfactored total loads).
for (i, m) in enumerate(sp.minimizers)
    s = I_symm(m...)
    sp.model.elements[:beam][i].section = Section(s.A, steel_ksi.E, steel_ksi.G,
                                                  s.Ix, s.Iy, s.J)
end
update_load_values_staged!(sp.model, sp, load_case=:all)
sp.load_dictionary = get_load_dictionary_by_id(sp.model)
Asap.solve!(sp.model, reprocess=true)
beam_ids = [get_element_id(be) for be in sp.model.elements[:beam]]
δ_bare_total_1 = let disp = ElementDisplacements(sp.model.elements[:beam][1],
                                                 sp.load_dictionary[beam_ids[1]],
                                                 resolution=200)
    maximum(abs.(disp.ulocal[2, :]))
end
println("  beam 1:  L = $(round(L1_in, digits=2)) in   L/360 = $(round(L1_in/360, digits=3)) in")
println("  beam 1:  bare Ix = $(round(bare_Ix1, digits=2)) in⁴")
println("  beam 1:  δ_bare(total unfactored) = $(round(δ_bare_total_1, digits=3)) in")
println("            ratio δ/(L/360) = $(round(δ_bare_total_1 / (L1_in/360), digits=2))")

println("\n── Step 3: staged breakdown (manual replay of postprocess_slab) ──")
# 3a: slab DL on bare Ix
update_load_values_staged!(sp.model, sp, load_case=:slab_dead)
sp.load_dictionary = get_load_dictionary_by_id(sp.model)
Asap.solve!(sp.model, reprocess=true)
δ_slab_dead_1 = let disp = ElementDisplacements(sp.model.elements[:beam][1],
                                                sp.load_dictionary[beam_ids[1]],
                                                resolution=200)
    maximum(abs.(disp.ulocal[2, :]))
end

# 3b: SDL + live on bare Ix (composite_action=false → no composite Ix boost)
update_load_values_staged!(sp.model, sp, load_case=:sdl_live)
sp.load_dictionary = get_load_dictionary_by_id(sp.model)
Asap.solve!(sp.model, reprocess=true)
δ_sdl_live_1 = let disp = ElementDisplacements(sp.model.elements[:beam][1],
                                               sp.load_dictionary[beam_ids[1]],
                                               resolution=200)
    maximum(abs.(disp.ulocal[2, :]))
end

# 3c: analytical self-weight on bare Ix
w_sw_1 = sec1.A * steel_ksi.ρ
δ_beam_dead_1 = 5 * w_sw_1 * L1_in^4 / (384 * steel_ksi.E * bare_Ix1)

δ_total_1 = δ_slab_dead_1 + δ_beam_dead_1 + δ_sdl_live_1
println("  δ_slab_dead  = $(round(δ_slab_dead_1, digits=3)) in")
println("  δ_beam_dead  = $(round(δ_beam_dead_1, digits=3)) in  (analytic)")
println("  δ_sdl+live   = $(round(δ_sdl_live_1,  digits=3)) in")
println("  δ_total      = $(round(δ_total_1,     digits=3)) in")
println("  L/360        = $(round(L1_in/360, digits=3)) in")
println("  L/240        = $(round(L1_in/240, digits=3)) in")

println("\n── Step 4: postprocess_slab result (same sections, composite_action=true) ──")
# Rebuild a fresh sizing_params in composite mode, inject NLP sections,
# call postprocess_slab, and see what it reports.
sp_c = SlabSizingParams(live_load=psf_to_ksi(50), superimposed_dead_load=psf_to_ksi(15),
                        slab_dead_load=0.0, live_factor=1.6, dead_factor=1.2,
                        beam_sizer=:continuous, nlp_solver=:MMA, max_depth=25.0,
                        beam_units=:in, serviceability_lim=360, collinear=false,
                        minimum_continuous=true, n_max_sections=0,
                        composite_action=true)
conv = convert_to_m[slab.slab_units] * 1 / convert_to_m[sp_c.beam_units]
sp_c.area         = slab.area * conv^2
sp_c.slab_depth_in = maximum(slab.slab_depths) * conv   # ← critical!
sp_c.i_perimeter  = Set(slab.i_perimeter)
sp_c.max_bay_span = maximum(slab.max_spans) * conv
sp_c.model        = get_scaled_model(slab, sp_c, conv)
sp_c.load_dictionary = get_load_dictionary_by_id(sp_c.model)
sp_c.minimizers = deepcopy(nlp_mins)
sp_c.ids        = deepcopy(nlp_ids)
sp_c.minimums   = zeros(Float64, n_beams)

res_c = postprocess_slab(slab, sp_c; check_collinear=false)

println("  n_L360_fail = $(res_c.n_L360_fail) / $n_beams")
println("  n_L240_fail = $(res_c.n_L240_fail) / $n_beams")
println("  max δ_total = $(round(res_c.max_δ_total, digits=3)) in")
println("  beam 1 staged (from postprocess_slab):")
println("    δ_slab_dead = $(round(res_c.δ_slab_dead[1], digits=3)) in")
println("    δ_beam_dead = $(round(res_c.δ_beam_dead[1], digits=3)) in")
println("    δ_sdl       = $(round(res_c.δ_sdl[1],      digits=3)) in")
println("    δ_live      = $(round(res_c.δ_live[1],     digits=3)) in")
println("    δ_total     = $(round(res_c.δ_total[1],    digits=3)) in")
