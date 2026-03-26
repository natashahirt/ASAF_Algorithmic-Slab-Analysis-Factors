# Physical Constants
const GRAVITY = 9.80665 # N/kg

# Material Densities (kg/m³)
const ρ_CONCRETE = 2400.0
const ρ_STEEL = 7850.0
const ρ_REBAR = 7850.0

# Material Densities (kip/in³)
const ρ_CONCRETE_KIPIN3 = kgm3_to_kipin3(ρ_CONCRETE)
const ρ_STEEL_KIPIN3 = kgm3_to_kipin3(ρ_STEEL)
const ρ_REBAR_KIPIN3 = kgm3_to_kipin3(ρ_REBAR)

# Embodied Carbon Coefficients
const ECC_STEEL = 1.22
const ECC_CONCRETE = 0.152 # from CLF
const ECC_REBAR = 0.854
# SFRM (spray-applied fireproofing), kgCO₂e/kg steel-equivalent mass basis for norm_mass_fireproofing
# Blaze-Shield II HS EPD (Isolatek / SM Transparency), declared unit 1000 kg, GWP IPCC total A1–A3 ≈ 1156.5 kgCO₂e
const ECC_FIREPROOFING = 1156.5 / 1000.0

# Big M
const BIG_M = 1e9