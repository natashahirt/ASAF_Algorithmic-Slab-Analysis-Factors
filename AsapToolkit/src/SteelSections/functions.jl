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
        error("Section '$name' (cleaned: '$key') was not found in AISC database.")
    end
end

# Geometry from the database is in mm / mm² / mm⁴ (etc.); Asap expects SI (m, m², m⁴).
const _MM2_TO_M2 = 1e-6
const _MM4_TO_M4 = 1e-12

"""
    toASAPframe(section::TorsionAllowed, E, G; ρ = 7850.0)

Build an `Asap.Section` from tabulated metric properties. `E` and `G` are in Pa; `ρ` is kg/m³.
Section areas and inertias are converted from mm units to SI.
"""
function toASAPframe(section::TorsionAllowed, E::Real, G::Real; ρ::Real = 7850.0)
    A = section.A * _MM2_TO_M2
    Ix = section.Ix * _MM4_TO_M4
    Iy = section.Iy * _MM4_TO_M4
    J = section.J * _MM4_TO_M4
    return Asap.Section(A, E, G, Ix, Iy, J, ρ)
end

function toASAPframe(name::String, E::Real, G::Real; ρ::Real = 7850.0)
    return toASAPframe(get_section(name), E, G; ρ = ρ)
end

function toASAPframe(name::String; E::Real = 200e9, G::Real = 77e9, ρ::Real = 7850.0)
    return toASAPframe(name, E, G; ρ = ρ)
end

"""
    toASAPtruss(section::AbstractSection, E; ρ = 7850.0)

`E` in Pa; `ρ` in kg/m³. Area is converted from mm² to m².
"""
function toASAPtruss(section::AbstractSection, E::Real; ρ::Real = 7850.0)
    A = section.A * _MM2_TO_M2
    return Asap.TrussSection(A, E, ρ)
end

function toASAPtruss(name::String, E::Real; ρ::Real = 7850.0)
    return toASAPtruss(get_section(name), E; ρ = ρ)
end
