"""
    get_deflection_constraint_nov2024(beam_params, beam_length, params) -> Function

Byte-faithful port of
`TributaryAreas @ a8dd2b98:slab_analysis/size_beams_utils.jl::get_deflection_constraint`,
the exact deflection constraint that governed the Nov 24 2024 sweep.

The modern `get_deflection_constraint` in `VariableBeamOptimizer/sizer/size_beams_utils.jl`
diverges from the Nov 2024 form in three ways that matter for reproduction:

  1. It calls `get_element_deflection(..., dead_load=0.0)` and then adds an
     analytical self-weight term `δ_sw = 5·w·L⁴/(384·E·Ix)` outside the FE.
     For a prismatic simply-supported beam with symmetric point loads the
     superposition is exact, but the location of `maximum(abs.(δ_local))`
     and the closed-form `δ_sw` midspan both live at the same node, so
     the sum is equivalent. For an asymmetric loading (non-centred point
     loads, different tributary widths left/right), the two maxima can
     shift to different FE nodes and the analytical addition overshoots.
  2. It uses `maximum(abs.(δ_local))` instead of `abs(minimum(δ_local))`.
     Equivalent for downward-only deflection; can diverge if any node
     rebounds upward (negative δ with large absolute value elsewhere).
  3. It applies `/ params.deflection_reduction_factor`; Nov 2024 had no
     such knob.

This port keeps the Nov 2024 signature and semantics exactly:

  - One deflection call: `get_element_deflection(beam_params, vars, material=steel_ksi)`
    with `dead_load` defaulted to `material.ρ`, so beam self-weight is
    included inside the FE integration.
  - `abs(minimum(δ_local))` against `L / params.serviceability_lim`.
  - No `deflection_reduction_factor` scaling.

Used by `process_discrete_beams_nov2024` and `process_continuous_beams_nov2024`.
"""
function get_deflection_constraint_nov2024(beam_params::FrameOptParams,
                                           beam_length::Real,
                                           params::SlabSizingParams)

    function single_section_deflection(vars::Vector)
        δ_local = Zygote.ignore() do
            get_element_deflection(beam_params, vars, material=steel_ksi)
        end

        δ_max = beam_length / params.serviceability_lim
        δ_min = minimum(δ_local)
        δ     = abs(δ_min)
        return [δ - δ_max]
    end

    return single_section_deflection
end
