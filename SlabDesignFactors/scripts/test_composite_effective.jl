"""
Unit tests for `get_I_composite_effective` — spatially-varying composite stiffness
with strip aggregation and perimeter beam support.

Run standalone:  include("SlabDesignFactors/scripts/test_composite_effective.jl")
Requires only TributaryAreas/composite/composite.jl (no model or Asap dependency).
"""

include("../../TributaryAreas/composite/composite.jl")

# ── Section & slab parameters ────────────────────────────────────────────────
h, w, tw, tf = 12.0, 6.0, 0.25, 0.4   # W-shape dimensions (in)
t_slab = 5.0                            # slab thickness (in)
E_s = 29000.0                           # steel modulus (ksi)
E_c = 57.0 * sqrt(4000.0)              # concrete modulus (ksi)
L_beam = 300.0                          # beam span (in)
uniform_width = 80.0                    # constant tributary width (in)

n_tests_passed = 0
n_tests_total = 0

# ── Test 1: Uniform width reproduces single-call result ──────────────────────
n_tests_total += 1
b_eff = get_b_eff(L_beam, uniform_width)
I_old = get_I_composite(h, w, tw, tf, t_slab, b_eff, E_s, E_c)

positions = collect(range(0.05, 0.95, length=50))
widths = fill(uniform_width, 50)
I_new = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                   L_beam, positions, widths)

if isapprox(I_old, I_new, rtol=1e-3)
    n_tests_passed += 1
    println("✓ Test 1 PASSED: uniform width — I_old=$(round(I_old, digits=2)), I_new=$(round(I_new, digits=2))")
else
    println("✗ Test 1 FAILED: uniform width — I_old=$(round(I_old, digits=2)), I_new=$(round(I_new, digits=2))")
end

# ── Test 2: Narrow at midspan → I_eff < I_avg ────────────────────────────────
n_tests_total += 1
widths_narrow_mid = [100.0 - 60.0 * (1 - 4*(p - 0.5)^2) for p in positions]

I_eff = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                   L_beam, positions, widths_narrow_mid)
avg_width = sum(widths_narrow_mid) / length(widths_narrow_mid)
I_avg = get_I_composite(h, w, tw, tf, t_slab, get_b_eff(L_beam, avg_width), E_s, E_c)

if I_eff < I_avg
    n_tests_passed += 1
    println("✓ Test 2 PASSED: narrow midspan — I_eff=$(round(I_eff, digits=2)) < I_avg=$(round(I_avg, digits=2))")
else
    println("✗ Test 2 FAILED: narrow midspan — I_eff=$(round(I_eff, digits=2)), I_avg=$(round(I_avg, digits=2))")
end

# ── Test 3: Wide at midspan → I_eff > I_avg ──────────────────────────────────
n_tests_total += 1
widths_wide_mid = [10.0 + 40.0 * (1 - 4*(p - 0.5)^2) for p in positions]

I_eff2 = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                    L_beam, positions, widths_wide_mid)
avg_width2 = sum(widths_wide_mid) / length(widths_wide_mid)
I_avg2 = get_I_composite(h, w, tw, tf, t_slab, get_b_eff(L_beam, avg_width2), E_s, E_c)

if I_eff2 > I_avg2
    n_tests_passed += 1
    println("✓ Test 3 PASSED: wide midspan — I_eff=$(round(I_eff2, digits=2)) > I_avg=$(round(I_avg2, digits=2))")
else
    println("✗ Test 3 FAILED: wide midspan — I_eff=$(round(I_eff2, digits=2)), I_avg=$(round(I_avg2, digits=2))")
end

# ── Test 4: Edge cases — single point and empty ──────────────────────────────
n_tests_total += 1
I_single = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                      L_beam, [0.5], [80.0])
I_bare = _Ix_I_symm(h, w, tw, tf)
I_empty = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                     L_beam, Float64[], Float64[])

if I_single > 0 && isapprox(I_empty, I_bare)
    n_tests_passed += 1
    println("✓ Test 4 PASSED: edge cases — I_single=$(round(I_single, digits=2)), I_empty=$(round(I_empty, digits=2)) == I_bare=$(round(I_bare, digits=2))")
else
    println("✗ Test 4 FAILED: edge cases — I_single=$(round(I_single, digits=2)), I_empty=$(round(I_empty, digits=2)), I_bare=$(round(I_bare, digits=2))")
end

# ── Test 5: Strip aggregation — two one-sided widths sum correctly ────────────
# Simulate an interior beam with identical strips from two adjacent cycles.
# Each side contributes half the total tributary width. The aggregated result
# should match a single call with the full width.
n_tests_total += 1
half_width = 40.0
full_width = 80.0

positions_side1 = collect(range(0.05, 0.95, length=50))
positions_side2 = collect(range(0.05, 0.95, length=50))
widths_side1 = fill(half_width, 50)
widths_side2 = fill(half_width, 50)

I_aggregated = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                          L_beam,
                                          vcat(positions_side1, positions_side2),
                                          vcat(widths_side1, widths_side2))

I_full = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                    L_beam, positions_side1, fill(full_width, 50))

if isapprox(I_aggregated, I_full, rtol=1e-3)
    n_tests_passed += 1
    println("✓ Test 5 PASSED: strip aggregation — I_agg=$(round(I_aggregated, digits=2)) ≈ I_full=$(round(I_full, digits=2))")
else
    println("✗ Test 5 FAILED: strip aggregation — I_agg=$(round(I_aggregated, digits=2)), I_full=$(round(I_full, digits=2))")
end

# ── Test 6: Perimeter beam gives lower stiffness than interior ────────────────
# Same tributary width, but perimeter cap is L/8 vs interior L/4.
# With trib_width=60 and L=300: interior b_eff = min(75,60) = 60,
# perimeter b_eff = min(37.5,60) = 37.5  →  I_perim < I_interior.
n_tests_total += 1
test_width = 60.0
positions_t6 = collect(range(0.05, 0.95, length=50))
widths_t6 = fill(test_width, 50)

I_interior = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                        L_beam, positions_t6, widths_t6;
                                        is_perimeter=false)
I_perimeter = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                         L_beam, positions_t6, widths_t6;
                                         is_perimeter=true)

if I_perimeter < I_interior
    n_tests_passed += 1
    println("✓ Test 6 PASSED: perimeter < interior — I_perim=$(round(I_perimeter, digits=2)) < I_int=$(round(I_interior, digits=2))")
else
    println("✗ Test 6 FAILED: perimeter < interior — I_perim=$(round(I_perimeter, digits=2)), I_int=$(round(I_interior, digits=2))")
end

# ── Test 7: Perimeter b_eff matches manual calculation ────────────────────────
# For uniform-width perimeter beam: b_eff = min(L/8, trib_width).
# With L=300, trib_width=50: b_eff = min(37.5, 50) = 37.5.
n_tests_total += 1
perim_width = 50.0
b_eff_expected = get_b_eff(L_beam, perim_width; is_perimeter=true)
I_expected_perim = get_I_composite(h, w, tw, tf, t_slab, b_eff_expected, E_s, E_c)

I_perim_effective = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                               L_beam,
                                               collect(range(0.05, 0.95, length=50)),
                                               fill(perim_width, 50);
                                               is_perimeter=true)

if isapprox(I_expected_perim, I_perim_effective, rtol=1e-3)
    n_tests_passed += 1
    println("✓ Test 7 PASSED: perimeter b_eff — I_expected=$(round(I_expected_perim, digits=2)), I_eff=$(round(I_perim_effective, digits=2))")
else
    println("✗ Test 7 FAILED: perimeter b_eff — I_expected=$(round(I_expected_perim, digits=2)), I_eff=$(round(I_perim_effective, digits=2))")
end

# ── Test 8: Asymmetric aggregation — unequal widths from two sides ────────────
# Interior beam where one cycle contributes 30 in and the other 60 in at each
# station. The aggregated total (90 in) should be used, not the average or max.
n_tests_total += 1
positions_asym = collect(range(0.05, 0.95, length=50))
widths_sideA = fill(30.0, 50)
widths_sideB = fill(60.0, 50)

I_asym = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                    L_beam,
                                    vcat(positions_asym, positions_asym),
                                    vcat(widths_sideA, widths_sideB))

I_total_90 = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                        L_beam, positions_asym, fill(90.0, 50))

if isapprox(I_asym, I_total_90, rtol=1e-3)
    n_tests_passed += 1
    println("✓ Test 8 PASSED: asymmetric aggregation — I_asym=$(round(I_asym, digits=2)) ≈ I_90=$(round(I_total_90, digits=2))")
else
    println("✗ Test 8 FAILED: asymmetric aggregation — I_asym=$(round(I_asym, digits=2)), I_90=$(round(I_total_90, digits=2))")
end

# ── Test 9: Interleaved strips from two sides (unequal counts) ────────────────
# One cycle places 20 strips, the other 30 strips — different positions, so most
# don't merge. Result should still be >= bare steel and <= the capped-uniform case.
n_tests_total += 1
pos_side1 = collect(range(0.05, 0.95, length=20))
pos_side2 = collect(range(0.04, 0.96, length=30))
w_side1 = fill(35.0, 20)
w_side2 = fill(45.0, 30)

I_interleaved = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                           L_beam,
                                           vcat(pos_side1, pos_side2),
                                           vcat(w_side1, w_side2))

I_capped = get_I_composite(h, w, tw, tf, t_slab, get_b_eff(L_beam, 200.0), E_s, E_c)

if I_interleaved > I_bare && I_interleaved <= I_capped
    n_tests_passed += 1
    println("✓ Test 9 PASSED: interleaved strips — I_bare=$(round(I_bare, digits=2)) < I_eff=$(round(I_interleaved, digits=2)) ≤ I_cap=$(round(I_capped, digits=2))")
else
    println("✗ Test 9 FAILED: interleaved strips — I_bare=$(round(I_bare, digits=2)), I_eff=$(round(I_interleaved, digits=2)), I_cap=$(round(I_capped, digits=2))")
end

# ── Test 10: Very sparse strips (2–3 points) ─────────────────────────────────
# Short beam or coarse spacing — make sure dx logic doesn't blow up.
n_tests_total += 1
I_2pt = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                   L_beam, [0.25, 0.75], [60.0, 60.0])
I_3pt = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                   L_beam, [0.2, 0.5, 0.8], [50.0, 70.0, 50.0])

if I_2pt > I_bare && I_3pt > I_bare
    n_tests_passed += 1
    println("✓ Test 10 PASSED: sparse strips — I_2pt=$(round(I_2pt, digits=2)), I_3pt=$(round(I_3pt, digits=2)) > I_bare=$(round(I_bare, digits=2))")
else
    println("✗ Test 10 FAILED: sparse strips — I_2pt=$(round(I_2pt, digits=2)), I_3pt=$(round(I_3pt, digits=2)), I_bare=$(round(I_bare, digits=2))")
end

# ── Test 11: All widths exceed L/4 cap → saturates to capped-uniform ─────────
# When every strip width >> L/4, get_b_eff clips them all to L/4. The result
# should equal the uniform case at the cap width.
n_tests_total += 1
huge_width = 500.0  # >> L/4 = 75
I_huge = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                    L_beam, positions, fill(huge_width, 50))
b_eff_cap = get_b_eff(L_beam, huge_width)  # = L/4 = 75
I_at_cap = get_I_composite(h, w, tw, tf, t_slab, b_eff_cap, E_s, E_c)

if isapprox(I_huge, I_at_cap, rtol=1e-3)
    n_tests_passed += 1
    println("✓ Test 11 PASSED: width saturation — I_huge=$(round(I_huge, digits=2)) ≈ I_cap=$(round(I_at_cap, digits=2))")
else
    println("✗ Test 11 FAILED: width saturation — I_huge=$(round(I_huge, digits=2)), I_cap=$(round(I_at_cap, digits=2))")
end

# ── Test 12: I_eff always >= bare steel Ix (sanity bound) ─────────────────────
# Composite action can only add stiffness. Check across several configurations.
n_tests_total += 1
configs_pass = true
for (pos_vec, w_vec, perim) in [
    (positions, fill(20.0, 50), false),
    (positions, fill(20.0, 50), true),
    ([0.5], [10.0], false),
    (collect(range(0.1, 0.9, length=5)), [15.0, 30.0, 60.0, 30.0, 15.0], false),
    (collect(range(0.1, 0.9, length=5)), [15.0, 30.0, 60.0, 30.0, 15.0], true),
]
    I_test = get_I_composite_effective(h, w, tw, tf, t_slab, E_s, E_c,
                                        L_beam, pos_vec, w_vec; is_perimeter=perim)
    if I_test < I_bare - 1e-6
        global configs_pass = false
        println("  ✗ I_eff=$(round(I_test, digits=2)) < I_bare=$(round(I_bare, digits=2)) for perim=$perim, widths=$w_vec")
    end
end

if configs_pass
    n_tests_passed += 1
    println("✓ Test 12 PASSED: I_eff >= I_bare for all configurations")
else
    println("✗ Test 12 FAILED: some configuration produced I_eff < I_bare")
end

# ── Summary ──────────────────────────────────────────────────────────────────
println("\n$(n_tests_passed)/$(n_tests_total) tests passed.")
if n_tests_passed == n_tests_total
    println("All tests passed.")
else
    error("$(n_tests_total - n_tests_passed) test(s) failed!")
end
