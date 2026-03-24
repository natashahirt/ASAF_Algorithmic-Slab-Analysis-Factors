#get all tabulated sections
allW() = [W(name) for name in names[Wrange]]
allC() = [C(name) for name in names[Crange]]
allL() = [L(name) for name in names[Lrange]]
allLL() = [LL(name) for name in names[LLrange]]
allWT() = [WT(name) for name in names[WTrange]]
allHSSRect() = [HSSRect(name) for name in names[HSSRectrange]]
allHSSRound() = [HSSRound(name) for name in names[HSSRoundrange]]

#get all names
Wnames = names[Wrange]
Cnames = names[Crange]
Lnames = names[Lrange]
LLnames = names[LLrange]
WTnames = names[WTrange]
HSSRectnames = names[HSSRectrange]
HSSRoundnames = names[HSSRoundrange]

# cache access to sections
const SECTION_CACHE = Dict{String, AbstractSection}()

# populate cache (fuzzy matching helper)
clean_name(n) = replace(uppercase(n), " " => "")

function populate_cache!()
    empty!(SECTION_CACHE)
    
    function cache_section!(s)
        SECTION_CACHE[clean_name(s.name)] = s
        SECTION_CACHE[clean_name(s.name_imperial)] = s
    end

    foreach(cache_section!, allW())
    foreach(cache_section!, allC())
    foreach(cache_section!, allL())
    foreach(cache_section!, allLL())
    foreach(cache_section!, allWT())
    foreach(cache_section!, allHSSRect())
    foreach(cache_section!, allHSSRound())
    
    return nothing
end

"""
get a steel section by name (metric or imperial)
is case-insensitive and space-insensitive.
"""
function get_section(name::String)
    # auto-populates the cache on first call.
    if isempty(SECTION_CACHE)
        populate_cache!()
    end
    
    key = clean_name(name)
    if haskey(SECTION_CACHE, key)
        return SECTION_CACHE[key]
    else
        error("Section '$name' (cleaned: '$key') not found in AISC database.")
    end
end

# convert units using Unitful - returns Unitful Section
function toASAPframe(section::TorsionAllowed, E::Quantity, G::Quantity; unit = u"mm", ρ = 7850.0u"kg/m^3")
    # Section properties already have Unitful units (mm², mm⁴, etc.)
    # Convert directly to SI - no need to multiply by unit^n
    A = uconvert(u"m^2", section.A)
    E_si = uconvert(u"Pa", E)
    G_si = uconvert(u"Pa", G)
    Ix = uconvert(u"m^4", section.Ix)
    Iy = uconvert(u"m^4", section.Iy)
    J = uconvert(u"m^4", section.J)
    ρ_si = uconvert(u"kg/m^3", ρ)
    
    return Asap.Section(A, E_si, G_si, Ix, Iy, J, ρ_si)
end

# only using section name
function toASAPframe(name::String, E::Quantity, G::Quantity; unit = u"mm", ρ = 7850.0u"kg/m^3")
    return toASAPframe(get_section(name), E, G; unit = unit, ρ = ρ)
end

# including standard defaults
function toASAPframe(name::String; E = 200u"GPa", G = 77u"GPa", unit = u"mm", ρ = 7850.0u"kg/m^3")
    return toASAPframe(name, E, G; unit = unit, ρ = ρ)
end

# fallback for Real (promotes to Unitful)
function toASAPframe(section::TorsionAllowed, E::Real, G::Real; unit = u"mm", ρ::Real = 7850.0)
    return toASAPframe(section, E * u"Pa", G * u"Pa"; unit = unit, ρ = ρ * u"kg/m^3")
end

# convert units using Unitful - returns Unitful TrussSection
function toASAPtruss(section::AbstractSection, E::Quantity; unit = u"mm", ρ = 7850.0u"kg/m^3")
    # Section properties already have Unitful units - convert directly to SI
    A = uconvert(u"m^2", section.A)
    E_si = uconvert(u"Pa", E)
    ρ_si = uconvert(u"kg/m^3", ρ)
    
    return Asap.TrussSection(A, E_si, ρ_si)
end

# fallback for Real (promotes to Unitful)
function toASAPtruss(section::AbstractSection, E::Real; unit = u"mm", ρ::Real = 7850.0)
    return toASAPtruss(section, E * u"Pa"; unit = unit, ρ = ρ * u"kg/m^3")
end
