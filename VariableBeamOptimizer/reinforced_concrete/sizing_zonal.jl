mutable struct LongitudinalZone

    index_range::Vector{Float64} # tuple, two indices 
    x_range::Vector{Float64} # tuple, two x values 
    Mn::Vector{Float64} # moment in that zone
    sign::Int64
    max_Mn::Float64
    As::Float64
    As′::Float64
    actual_As::Float64
    actual_As′::Float64
    rebar_top::Vector{Vector{Int64}}
    rebar_bottom::Vector{Vector{Int64}}

    function LongitudinalZone(index_range, x_range, Mn)
        max_Mn = maximum(abs.(Mn))
        sign = Mn[1] >= 0 ? 1 : -1
        As = As′ = 0
        actual_As = actual_As′ = 0.0
        rebar_top = Vector{Vector{Int64}}()      # []
        rebar_bottom = Vector{Vector{Int64}}()   # []
        new(index_range, x_range, Mn, sign, max_Mn, As, As′, actual_As, actual_As′, rebar_top, rebar_bottom)
    end

end

"""
structure to store stirrup information for each stirrup zone (i.e. changes in shear magnitude)
"""
mutable struct StirrupZone
    spacing::Float64               # stirrup spacing in inches
    index_range::UnitRange{Int}    # index range in x_positions
    x_range::Tuple{Float64, Float64} # physical x extent of the zone
    stirrups::Vector{Float64}     # x positions of valid stirrups
    n_stirrups::Int64              # number of stirrups

    function StirrupZone(spacing, index_range, x_range, stirrups::Vector{Float64}=Float64[])
        n_stirrups = length(stirrups)
        new(spacing, index_range, x_range, stirrups, n_stirrups)
    end
end

"""
longitudinal zones are selected based on when the moment changes sign
since the larger rebar needs to be resisting tension
"""
function get_longitudinal_zones(Mn::Vector{Float64}, x::Vector{Float64})
    zone_start = 1
    zones = LongitudinalZone[]

    for i in 2:lastindex(Mn)
        if i==lastindex(Mn) || (Mn[i-1]) * Mn[i] < 0
            if i == lastindex(Mn)
                i+=1
            end
            push!(zones, LongitudinalZone([zone_start, i-1], 
                                        [x[zone_start], x[i-1]],
                                        Mn[zone_start:i-1]))
            zone_start = i
        end
    end
    return zones, nothing
end

function get_longitudinal_zones(beams::Vector{Vector{Float64}}, x::Vector{Float64})

    zone_start = 1
    zones = LongitudinalZone[]
    beam_ids = Int[]

    for (idx, Mn) in enumerate(beams)
        for i in 2:lastindex(Mn)
            if i==lastindex(Mn) || (Mn[i-1]) * Mn[i] < 0
                if i == lastindex(Mn)
                    i+=1
                end
                new_zone = LongitudinalZone(  [zone_start, i-1], 
                                                [x[zone_start], x[i-1]],
                                                Mn[zone_start:i-1])
                push!(zones, new_zone)
                push!(beam_ids, idx)
                zone_start = i
            end
        end
    end

    return zones, beam_ids
end

"""
Longitudinal as well as stirrup steel needs to be treated zone-by-zone
Longitudinal
    - distinguish between negative and positive moment regions
    - assume #3 rebar for stirrups
"""
function size_concrete_zonal(b, d, x::Vector{Float64}, beams::Vector{Float64}, Vn::Vector{Float64}; material::ReinforcedConcrete=ReinforcedConcrete(5,60), ϕb::Float64 = 0.9, ϕv::Float64 = 0.75, d′=1., clear_side=1., clear_bottom=2.5, clear_top=1., plot::Bool=false)
    zones, beam_ids = get_longitudinal_zones(beams, x)

    for zone in zones
        Mu = ϕb * zone.max_Mn
        zone.As, zone.As′ = get_longitudinal_As_differentiable(b, d, d′, Mu, material)
    end

    total_As = select_rebar_gurobi(b, d, zones, material)

    return total_As

end