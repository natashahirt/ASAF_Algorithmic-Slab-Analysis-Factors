"""
    Metal(E::Real, G::Real, Fy::Real, Fu::Real, ρ::Real, ν::Real, units::Symbol)

Structural metal (not to be confused with Asap's Material)
If you use ksi, geometric units should be in inches. 
If you use MPa, geometric units should be in mm.

# Arguments
- E::Real : Young's modulus [ksi, MPa]
- G::Real : shear modulus [ksi, MPa]
- Fy::Real : yield strength [ksi, MPa]
- Fu::Real :  ultimate strength [ksi, MPa]
- ρ::Real : volumetric density [g/mm³, g/cm³, kip/in³]
- ν::Real : Poisson's ratio in elastic range [-,-]
- units::Symbol : what units do we need for the geometric stuff? [in, mm]
"""
struct Metal <: AbstractMaterial
    E::Real #young's modulus
    G::Real #shear modulus
    Fy::Real # yield strength
    Fu::Real # ultimate strength
    ρ::Real #density 
    ν::Real #poisson's ratio
    units::Symbol # what units are we using

    function Metal(E::Real, G::Real, Fy::Real, Fu::Real, ρ::Real, ν::Real, units::Symbol)

        return new(E, G, Fy, Fu, ρ, ν, units)
        
    end
    
end

const steel_mm = Metal(200e3, 79.3e3, 250, 482, 0.00785, 0.26, :mm) # N/mm²
const steel_ksi = Metal(29e3, 11.5e3, 50, 65, 0.2836 * 1e-3, 0.26, :in) # ksi, density in kip/in^3
const steel_cm = Metal(200e5, 79.3e5, 250e2, 482e2, 7.85, 0.26, :cm) # N/cm²
const steel_m = Metal(200e9, 79.3e9, 250e6, 482e6, 7850.0, 0.26, :m) # N/m²
const rebar_40_ksi = Metal(29e3, 11.5e3, 40, 70, 0.2836 * 1e-3, 0.26, :in) # ksi, density in kip/in^3
const rebar_60_ksi = Metal(29e3, 11.5e3, 60, 90, 0.2836 * 1e-3, 0.26, :in) # ksi, density in kip/in^3
const rebar_75_ksi = Metal(29e3, 11.5e3, 75, 100, 0.2836 * 1e-3, 0.26, :in) # ksi, density in kip/in^3

const material_dict = Dict(
    :in => steel_ksi,
    :m => steel_m,
    :mm => steel_mm,
    :cm => steel_cm
)