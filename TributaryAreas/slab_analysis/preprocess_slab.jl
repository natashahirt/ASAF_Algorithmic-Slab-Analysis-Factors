"""
    get_slab_self_weight(cycles, elements, slab_type; ρ_CONCRETE=2400.0, analysis=:cellular, precision=4, vector_1d=[1.0, 0.0])

Calculate the self-weight of a slab based on its type and analysis method.

# Arguments
- `cycles`: The cycles of the slab.
- `elements`: The elements of the slab.
- `slab_type`: The type of slab (:isotropic, :uniaxial, :orth_biaxial, :orthogonal).
- `analysis`: Analysis method (:cellular or :uniform).
- `precision`: Number of decimal places for rounding (default: 4).
- `vector_1d`: Direction vector for uniaxial slabs (default: [1.0, 0.0]).

# Returns
- `max_spans`: Maximum spans for each cycle.
- `slab_depths`: Calculated slab depths.
- `slab_loads`: Slab loads in kN/m².
"""
function _apply_slab_depth_minimum(self::SlabAnalysisParams, slab_depths::Vector{<:Real})::Vector{Float64}
    _floor = 0.001 / convert_to_m[self.slab_units]
    dmin = max(Float64(self.slab_depth_minimum), _floor)
    return Float64[max(Float64(d), dmin) for d in slab_depths]
end

function get_slab_depths(self::SlabAnalysisParams, cycles, elements; precision::Int64=4)
    max_spans = calculate_max_spans(self, cycles, elements)

    if self.slab_thickness > 0.0
        slab_depths = [self.slab_thickness for _ in cycles]
    else
        slab_depths = calculate_slab_depths(max_spans, self.slab_type)
    end
    slab_depths = _apply_slab_depth_minimum(self, slab_depths)

    if self.slab_sizer == :cellular
        print_results(max_spans, slab_depths, precision)
        return max_spans, slab_depths
    elseif self.slab_sizer == :uniform
        return final_get_slab_depths(cycles, max_spans, slab_depths, precision)
    end
end

function calculate_max_spans(self::SlabAnalysisParams, cycles, elements)
    n = length(cycles)
    max_spans = Vector{Float64}(undef, n)
    for (idx, cycle) in enumerate(cycles)
        startpoints, stiffnesses = get_start_and_end_coords(cycle, elements)
        max_span = if self.slab_type == :isotropic 
            get_max_diagonal_span_2d(self, startpoints, method=:orthogonal)
        elseif self.slab_type == :uniaxial
            cycle_edges = get_cycle_edges(startpoints)
            get_max_diagonal_span_1d(cycle_edges, self.vector_1d)
        else
            get_max_diagonal_span_2d(self, startpoints, method=:orthonormal, vector_1d=self.vector_1d)
        end
        validate_max_span(max_span, self.slab_type)
        max_spans[idx] = max_span
    end
    return max_spans
end

function validate_max_span(max_span, slab_type)
    max_allowed_span = slab_type == :uniaxial ? 12.5 : 13.6
    @assert max_span <= max_allowed_span "Your maximum span ($max_span m) exceeds the recommended maximum span for a $(slab_type) slab ($max_allowed_span m)"
end

function calculate_slab_depths(max_spans, slab_type)
    divisor = slab_type == :uniaxial ? 28 : 36
    return max_spans ./ divisor
end

function print_results(max_spans, slab_depths, precision)
    println("\n cycle maximum spans: $([round(max_span, digits=precision) for max_span in max_spans])")
    println(" slab thicknesses: $([round(slab_depth, digits=precision) for slab_depth in slab_depths])")
    println("\n")
end

function final_get_slab_depths(cycles, max_spans, slab_depths, precision)
    max_span = maximum(max_spans)
    max_calculated_depth = maximum(slab_depths)
    slab_depths = fill(max_calculated_depth, length(cycles))

    println("\n maximum span: $(round(max_span, digits=precision))")
    println(" slab thickness: $(round(max_calculated_depth, digits=precision))\n")

    return max_spans, slab_depths
end
