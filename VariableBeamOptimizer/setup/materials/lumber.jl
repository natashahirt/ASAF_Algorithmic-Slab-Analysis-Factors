"""
    Lumber(E::Real, G::Real, Fx::Real, Fxy::Real, ρ::Real, ν::Real, units::Symbol)

Structural lumber material properties.
Values should be provided in consistent units based on the units parameter.

# Arguments
- E::Real : Mean modulus of elasticity [N/mm², N/m²]
- G::Real : Shear modulus [N/mm², N/m²] 
- Fx::Real : Bending strength parallel to grain [N/mm², N/m²]
- Fxy::Real : Shear strength parallel to grain [N/mm², N/m²]
- ρ::Real : Mean density [kg/m³]
- ν::Real : Poisson's ratio [-]
- units::Symbol : Unit system for geometric properties [:mm, :m]
"""
struct Lumber <: AbstractMaterial
    E::Real # Mean modulus of elasticity
    G::Real # Shear modulus
    Fx::Real # Bending strength parallel to grain
    Fxy::Real # Shear strength parallel to grain
    ρ::Real # Mean density
    ν::Real # Poisson's ratio
    units::Symbol # Geometric unit system
    
    # Additional properties
    E_min::Real # Minimum modulus of elasticity
    ρ_005::Real # 5th percentile density
    ft_parallel::Real # Tension strength parallel to grain
    fc_parallel::Real # Compression strength parallel to grain  
    Fy::Real # Compression/tension strength perpendicular to grain

    function Lumber(E::Real, G::Real, Fx::Real, Fxy::Real, ρ::Real, ν::Real, units::Symbol,
                   E_min::Real, ρ_005::Real, ft_parallel::Real, fc_parallel::Real, Fy::Real)
        return new(E, G, Fx, Fxy, ρ, ν, units, E_min, ρ_005, ft_parallel, fc_parallel, Fy)
    end
end

# Common lumber grades based on provided strength data
const C16_mm = Lumber(
    8800.0, # E - Mean modulus of elasticity (N/mm²)
    360.0,  # G (N/mm²)
    5.3,    # Fx - Bending strength parallel to grain (N/mm²)
    0.67,   # Fxy - Shear strength parallel (N/mm²)
    370.0,  # ρ - Mean density (kg/m³)
    0.3,    # ν
    :mm,    # units
    5800.0, # E_min - Minimum modulus of elasticity (N/mm²)
    310.0,  # ρ_005 - 5th percentile density (kg/m³)
    3.2,    # ft_parallel - Tension strength parallel to grain (N/mm²)
    1.8,    # fc_parallel - Compression strength parallel to grain (N/mm²)
    2.2     # Fy - Compression strength perpendicular to grain (N/mm²)
) 

const C24_mm = Lumber(
    10800.0, # E - Mean modulus of elasticity (N/mm²)
    440.0,   # G (N/mm²)
    7.5,     # Fx - Bending strength parallel to grain (N/mm²)
    0.71,    # Fxy - Shear strength parallel (N/mm²)
    420.0,   # ρ - Mean density (kg/m³)
    0.3,     # ν
    :mm,     # units
    7200.0,  # E_min - Minimum modulus of elasticity (N/mm²)
    350.0,   # ρ_005 - 5th percentile density (kg/m³)
    4.5,     # ft_parallel - Tension strength parallel to grain (N/mm²)
    7.9,     # fc_parallel - Compression strength parallel to grain (N/mm²)
    2.4      # Fy - Compression strength perpendicular to grain (N/mm²)
)

const TR26_mm = Lumber(
    11000.0, # E - Mean modulus of elasticity (N/mm²)
    670.0,   # G (N/mm²)
    28.3,    # Fx - Bending strength parallel to grain (N/mm²)
    4.0,     # Fxy - Shear strength parallel (N/mm²)
    444.0,   # ρ - Mean density (kg/m³)
    0.3,     # ν
    :mm,     # units
    7400.0,  # E_min - Minimum modulus of elasticity (N/mm²)
    370.0,   # ρ_005 - 5th percentile density (kg/m³)
    17.6,    # ft_parallel - Tension strength parallel to grain (N/mm²)
    22.9,    # fc_parallel - Compression strength parallel to grain (N/mm²)
    2.6      # Fy - Compression strength perpendicular to grain (N/mm²)
)

const lumber_dict = Dict(
    :C16 => C16_mm,
    :C24 => C24_mm,
    :TR26 => TR26_mm
)