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

# Embodied Carbon Coefficients (kgCO₂e/kg, stages A1–A3, cradle-to-gate)
# Source: 2023 CLF North American Material Baselines (Waldman et al. 2023).
const ECC_STEEL    = 1.22   # Hot-rolled sections, fabricated (AISC 2021 IW-EPD)
const ECC_CONCRETE = 0.152  # Ready-mixed concrete, 5000 psi (34.5 MPa), US national average (365 kgCO₂e/m³ ÷ 2400 kg/m³)
const ECC_REBAR    = 0.854  # Rebar, fabricated (CRSI 2022 IW-EPD)

# Big M
const BIG_M = 1e9