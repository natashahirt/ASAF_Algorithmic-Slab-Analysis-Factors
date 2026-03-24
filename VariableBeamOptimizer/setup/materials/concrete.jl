mutable struct Concrete <: AbstractMaterial
    f_c::Float64 # compressive strength
    f_t::Float64 # tensile strength
    f_r::Float64 # modulus of rupture
    E_c::Float64 # elastic modulus
    ν::Float64 # poisson's ratio
    ρ::Float64 # density
    β₁::Float64 # compression block factor
    units::Symbol # units of measurement
end

const concrete_4000_psi = Concrete(4000, 0.12 * 4000, 7.5 * sqrt(4000), 57000 * sqrt(4000), 0.2, 150, 0.85, :psi)
const concrete_4_ksi = Concrete(4.0, 0.12 * 4.0, 7.5 * sqrt(4.0), 57000 * sqrt(4.0), 0.2, 150, 0.85, :ksi)
const concrete_30_mpa = Concrete(30.0, 0.12 * 30.0, 0.62 * sqrt(30.0), 4700 * sqrt(30.0), 0.2, 2400, 0.85, :mpa)

mutable struct ReinforcedConcrete <: AbstractMaterial

    fc′::Float64
    fy::Float64

    function ReinforcedConcrete(fc′, fy) 
        new(fc′, fy)
    end

end