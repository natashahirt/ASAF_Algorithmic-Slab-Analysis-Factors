"""
    to_ASAP_section(section::AbstractSection)

Turns variable section into an ASAP section (extracts section properties and puts them into the
ASAP struct format.)
Having material be custom-entered allows us to switch around units.
"""
function to_ASAP_section(section::AbstractSection; material::Union{Nothing,AbstractMaterial,Asap.Material}=nothing)

    if isnothing(material)
        material = section.material
    end

    return Section(section.A, material.E, material.G, section.Ix, section.Iy, section.J)

end