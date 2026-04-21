const etype2DOF = Dict(
    Element{Asap.FixedFixed} => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    Element{Asap.FreeFree} => [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    Element{Asap.FixedFree} => [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    Element{Asap.FreeFixed} => [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0],
    Element{Asap.Joist} => [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0]
)

"""End-release type parameter `T` for `Asap.Element{T}` (e.g. `Asap.FixedFixed`)."""
get_release_type(::Asap.Element{T}) where {T} = T

const release_type_to_symbol = Dict(
    Asap.FixedFixed => :fixedfixed,
    Asap.FreeFree => :freefree,
    Asap.FixedFree => :fixedfree,
    Asap.FreeFixed => :freefixed,
    Asap.Joist => :joist,
)

"""Symbol accepted by `Asap.Element(...; release=...)` for a given element."""
asap_release_symbol(el::Asap.Element) = release_type_to_symbol[get_release_type(el)]

const release2DOF = Dict(
    Asap.FixedFixed => etype2DOF[Element{Asap.FixedFixed}],
    Asap.FreeFree => etype2DOF[Element{Asap.FreeFree}],
    Asap.FixedFree => etype2DOF[Element{Asap.FixedFree}],
    Asap.FreeFixed => etype2DOF[Element{Asap.FreeFixed}],
    Asap.Joist => etype2DOF[Element{Asap.Joist}],
)

const MPointLoad = Dict(
    Asap.FixedFixed => MPoint_fixedfixed,
    Asap.FreeFree => MPoint_freefree,
    Asap.FixedFree => MPoint_fixedfree,
    Asap.FreeFixed => MPoint_freefixed,
    Asap.Joist => MPoint_freefree,
)

const VPointLoad = Dict(
    Asap.FixedFixed => VPoint_fixedfixed,
    Asap.FreeFree => VPoint_freefree,
    Asap.FixedFree => VPoint_fixedfree,
    Asap.FreeFixed => VPoint_freefixed,
    Asap.Joist => VPoint_freefree,
)

export planarDOFs
const planarDOFs = Dict(:X => [1, 7],
    :XY => [2, 6, 8, 12],
    :XZ => [3, 5, 9, 11])

"""
Collect all elements that are part of the same continuous member (ie collect shattered elements)
"""
function groupbyid(elements::Vector{<:Asap.AbstractElement})
    ichecked = Vector{Int64}()
    indices = Vector{Vector{Int64}}()

    elementids = getproperty.(elements, :elementID)

    for id in elementids
        
        in(id, ichecked) && continue

        igroup = findall(elementids .== id)

        push!(indices, igroup)
        push!(ichecked, id)
    end

    return indices
end