"""
    ConcreteMaterial

Concrete material scenario bundling density, modulus, compressive strength, and
embodied carbon coefficient.  All values are stored in **both** unit systems so
that downstream code never needs to convert at runtime.

# Fields – SI
- `name::String`           — human-readable label (e.g. "Normal-weight 4 ksi")
- `ρ_concrete::Float64`    — density [kg/m³]
- `f′c_MPa::Float64`       — compressive strength [MPa]
- `E_c_MPa::Float64`       — elastic modulus [MPa]
- `ecc_concrete::Float64`  — GWP [kgCO₂e / kg concrete]

# Fields – Imperial (pre-converted)
- `ρ_concrete_kipin3::Float64` — density [kip/in³]
- `f′c_ksi::Float64`           — compressive strength [ksi]
- `E_c_ksi::Float64`           — elastic modulus [ksi]
"""
struct ConcreteMaterial
    name::String

    # SI
    ρ_concrete::Float64       # kg/m³
    f′c_MPa::Float64          # MPa
    E_c_MPa::Float64          # MPa
    ecc_concrete::Float64     # kgCO₂e / kg

    # Imperial (pre-converted)
    ρ_concrete_kipin3::Float64 # kip/in³
    f′c_ksi::Float64           # ksi
    E_c_ksi::Float64           # ksi
end

const _MPA_PER_KSI = 6.89476

"""
    ConcreteMaterial(name; ρ_concrete, f′c_MPa, ecc_concrete)

Convenience constructor that derives E_c and imperial equivalents automatically.
E_c is computed as `4700 √f′c` (MPa) per ACI 318 (metric form).
"""
function ConcreteMaterial(name::String;
                          ρ_concrete::Real,
                          f′c_MPa::Real,
                          ecc_concrete::Real)
    E_c_MPa  = 4700.0 * sqrt(f′c_MPa)
    ρ_kipin3 = kgm3_to_kipin3(ρ_concrete)
    f′c_ksi  = f′c_MPa / _MPA_PER_KSI
    E_c_ksi  = E_c_MPa / _MPA_PER_KSI
    return ConcreteMaterial(name,
        Float64(ρ_concrete), Float64(f′c_MPa), E_c_MPa, Float64(ecc_concrete),
        ρ_kipin3, f′c_ksi, E_c_ksi)
end

# ── Predefined scenarios ─────────────────────────────────────────────────────

const CONCRETE_NORMAL_4KSI = ConcreteMaterial(
    "NW_4ksi";
    ρ_concrete  = 2400.0,       # kg/m³
    f′c_MPa     = 27.6,         # ≈ 4 ksi
    ecc_concrete = 0.152,       # kgCO₂e/kg  (CLF baseline)
)

const CONCRETE_NORMAL_5KSI = ConcreteMaterial(
    "NW_5ksi";
    ρ_concrete  = 2400.0,
    f′c_MPa     = 34.5,         # ≈ 5 ksi
    ecc_concrete = 0.163,       # slightly higher GWP for higher strength
)

const CONCRETE_LW_4KSI = ConcreteMaterial(
    "LW_4ksi";
    ρ_concrete  = 1840.0,       # kg/m³ (typical structural LWC)
    f′c_MPa     = 27.6,
    ecc_concrete = 0.208,       # higher GWP due to expanded-aggregate processing
)

const CONCRETE_LW_3KSI = ConcreteMaterial(
    "LW_3ksi";
    ρ_concrete  = 1760.0,       # kg/m³
    f′c_MPa     = 20.7,         # ≈ 3 ksi
    ecc_concrete = 0.195,
)

"""
All predefined concrete material scenarios, keyed by short name.
"""
const CONCRETE_SCENARIOS = Dict{String, ConcreteMaterial}(
    "NW_4ksi" => CONCRETE_NORMAL_4KSI,
    "NW_5ksi" => CONCRETE_NORMAL_5KSI,
    "LW_4ksi" => CONCRETE_LW_4KSI,
    "LW_3ksi" => CONCRETE_LW_3KSI,
)

"""
Default concrete material — matches the legacy global constants.
"""
const DEFAULT_CONCRETE = CONCRETE_NORMAL_4KSI
