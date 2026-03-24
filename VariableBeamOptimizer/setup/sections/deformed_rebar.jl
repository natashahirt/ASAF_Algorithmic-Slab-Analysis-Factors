mutable struct DeformedRebar
    size::Int64
    material::Metal
    diameter::Float64
    weight::Float64
    area::Float64
    units::Symbol
end

const rebar_3 = DeformedRebar(3, rebar_60_ksi, 0.375, 0.376, 0.11, :in)
const rebar_4 = DeformedRebar(4, rebar_60_ksi, 0.500, 0.668, 0.20, :in)
const rebar_5 = DeformedRebar(5, rebar_60_ksi, 0.625, 1.043, 0.31, :in)
const rebar_6 = DeformedRebar(6, rebar_60_ksi, 0.750, 1.502, 0.44, :in)
const rebar_7 = DeformedRebar(7, rebar_60_ksi, 0.875, 2.044, 0.60, :in)
const rebar_8 = DeformedRebar(8, rebar_60_ksi, 1.000, 2.670, 0.79, :in)
const rebar_9 = DeformedRebar(9, rebar_60_ksi, 1.128, 3.400, 1.00, :in)
const rebar_10 = DeformedRebar(10, rebar_60_ksi, 1.270, 4.303, 1.27, :in)
const rebar_11 = DeformedRebar(11, rebar_60_ksi, 1.410, 5.313, 1.56, :in)
const rebar_14 = DeformedRebar(14, rebar_60_ksi, 1.693, 7.650, 2.25, :in)
const rebar_18 = DeformedRebar(18, rebar_60_ksi, 2.257, 13.600, 4.00, :in)

const REBAR_TYPES = Dict(
            3 => rebar_3,
            4 => rebar_4,
            5 => rebar_5,
            6 => rebar_6,
            7 => rebar_7,
            8 => rebar_8,
            9 => rebar_9,
            10 => rebar_10,
            11 => rebar_11,
            14 => rebar_14,
            18 => rebar_18
)