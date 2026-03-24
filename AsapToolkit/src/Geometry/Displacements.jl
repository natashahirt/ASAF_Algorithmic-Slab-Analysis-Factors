"""
displacement function for the transverse translation and in-plane rotation for a GIVEN PLANE IN THE LCS OF AN ELEMENT:

u_xy= [ v₁
        θ₁
        v₂
        θ₂ ]

IE for the local XY plane:

v₁, v₂ are the start and end displacements in the local Y direction
θ₁, θ₂ are the start and end rotations in the **local Z** direction (ie rotation in plane of local XY)

Gives:

v_y(x) = N × u_xy (translational displacement in local Y at point x)

"""
function N(x::Float64, L::Float64)
    n1 = 1 - 3(x/L)^2 + 2(x/L)^3
    n2 = x * (1 - x/L)^2
    n3 = 3(x/L)^2 - 2(x/L)^3
    n4 = x^2/L * (-1 + x/L)

    return [n1 n2 n3 n4]
end

"""
Axial displacement function: linear interpolation between start and end displacements
"""
function Naxial(x::Float64, L::Float64)
    n1 = 1 - x/L
    n2 = x / L

    return [n1 n2]
end

"""
    displacements(element::Element; n::Integer = 20)

Get the [3 × n] matrix where each column represents the local [x,y,z] displacement of the element from end forces
"""
function unodal(element::Element; n::Integer = 20)

    # Strip units from displacements for matrix multiplication (R and K are Float64)
    # DOFs 1-3 are translations (m), DOFs 4-6 are rotations (rad)
    ustart = Asap.to_displacement_vec(element.nodeStart.displacement)
    uend = Asap.to_displacement_vec(element.nodeEnd.displacement)
    
    # base properties
    ulocal = element.R * [ustart; uend] .* etype2DOF[typeof(element)]
    L = ustrip(u"m", element.length)

    # extracting relevant nodal DOFs
    uX = ulocal[[1, 7]]
    uY = ulocal[[2, 6, 8, 12]]
    uZ = ulocal[[3, 5, 9, 11]] .* [1, -1, 1, -1]

    # discretizing length of element
    xrange = range(0, L, n)

    nA = vcat(Naxial.(xrange, L)...)
    nT = vcat(N.(xrange, L)...)

    dx = nA * uX
    dy = nT * uY
    dz = nT * uZ

    # [dx' ; dy' ; dz']
    dx, dy, dz
end

"""
Accumlate the internal forces cause by a line load to an element
"""
function accumulatedisp!(
    load::LineLoad, 
    xvals::Vector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    # Strip units from length and section properties (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", load.element.length))
    E = ustrip(u"Pa", uconvert(u"Pa", load.element.section.E))
    Istrong = ustrip(u"m^4", uconvert(u"m^4", load.element.section.Ix))
    Iweak = ustrip(u"m^4", uconvert(u"m^4", load.element.section.Iy))

    # distributed load magnitudes in LCS (strip units from load.value)
    value_stripped = [ustrip(u"N/m", uconvert(u"N/m", v)) for v in load.value]
    wx, wy, wz = (R * value_stripped) .* [1, -1, -1]

    # Use dispatched DLine function based on element type
    Dy .-= DLine.(Ref(load.element), wy, L, xvals, E, Istrong)
    Dz .-= DLine.(Ref(load.element), wz, L, xvals, E, Iweak)
end

"""
Accumlate the internal forces cause by a point load to an element
"""
function accumulatedisp!(
    load::PointLoad, 
    xvals::Vector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    # Strip units from length and section properties (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", load.element.length))
    E = ustrip(u"Pa", uconvert(u"Pa", load.element.section.E))
    Istrong = ustrip(u"m^4", uconvert(u"m^4", load.element.section.Ix))
    Iweak = ustrip(u"m^4", uconvert(u"m^4", load.element.section.Iy))
    frac = load.position

    # distributed load magnitudes in LCS (strip units from load.value)
    value_stripped = [ustrip(u"N", uconvert(u"N", v)) for v in load.value]
    px, py, pz = (R * value_stripped) .* [1, -1, -1]

    # Use dispatched DPoint function based on element type
    Dy .-= DPoint.(Ref(load.element), py, L, xvals, frac, E, Istrong)
    Dz .-= DPoint.(Ref(load.element), pz, L, xvals, frac, E, Iweak)
end

"""
Accumulate the displacement cause by a GravityLoad (self-weight) to an element.
GravityLoad applies a distributed load equal to ρ * A * g in the global -Z direction.
"""
function accumulatedisp!(
    load::Asap.GravityLoad, 
    xvals::Vector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64})

    element = load.element
    R = element.R[1:3, 1:3]
    # Strip units from length and section properties
    L = ustrip(u"m", element.length)
    E = ustrip(u"Pa", element.section.E)
    Istrong = ustrip(u"m^4", element.section.Ix)
    Iweak = ustrip(u"m^4", element.section.Iy)

    # Self-weight per unit length: w = ρ * A * g (in N/m)
    ρ = ustrip(u"kg/m^3", element.section.ρ)
    A = ustrip(u"m^2", element.section.A)
    g = ustrip(u"m/s^2", load.factor)
    w_mag = ρ * A * g  # N/m

    # Global load vector: self-weight acts in -Z direction
    w_global = [0.0, 0.0, -w_mag]

    # Transform to local coordinate system
    wx, wy, wz = (R * w_global) .* [1, -1, -1]

    # Use dispatched DLine function based on element type
    Dy .-= DLine.(Ref(element), wy, L, xvals, E, Istrong)
    Dz .-= DLine.(Ref(element), wz, L, xvals, E, Iweak)
end

"""
    ulocal(element::Element, model::Model; resolution = 20)

Get the [3 × resolution] matrix of xyz displacements in LCS
"""
function ulocal(element::Element, model::Model; resolution = 20)
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", element.length))

    xinc = collect(range(0, L, resolution))

    D = unodal(element; n = resolution)

    element_loads = get_elemental_loads(model)

    for load in element_loads[element.elementID]
        accumulatedisp!(load, xinc, D[2,:], D[3,:])
    end

    return D
end

"""
    uglobal(element::Element, model::Model; resolution = 20)

Get the [3 × resolution] matrix of xyz displacements in GCS
"""
function uglobal(element::Element, model::Model; resolution = 20)
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", element.length))

    xinc = collect(range(0, L, resolution))

    D = unodal(element; n = resolution)

    element_loads = get_elemental_loads(model)

    for load in element_loads[element.elementID]
        accumulatedisp!(load, xinc, D[2,:], D[3,:])
    end

    return hcat([sum(Δ .* element.LCS) for Δ in eachcol(D)]...)
end

struct ElementDisplacements
    element::Element
    resolution::Integer
    x::Vector{Float64}
    ulocal::Matrix{Float64}
    uglobal::Matrix{Float64}
    basepositions::Matrix{Float64}
end

"""
    ElementDisplacements(element::Element, model::Model; resolution = 20)

Get the local/global displacements of an element
"""
function ElementDisplacements(element::Asap.AbstractElement, model::Asap.Model; resolution = 20)
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", element.length))

    xinc = collect(range(0, L, resolution))

    Dx, Dy, Dz = unodal(element; n = resolution)

    element_loads = get_elemental_loads(model)

    for load in element_loads[element.elementID]
        accumulatedisp!(load, xinc, Dy, Dz)
    end

    D = [Dx'; Dy'; Dz']

    Dglobal = hcat([sum(Δ .* element.LCS) for Δ in eachcol(D)]...)

    # Strip units from position (convert Quantity → Float64) for basepoints calculation
    pos_start = [ustrip(u"m", uconvert(u"m", p)) for p in element.nodeStart.position]
    basepoints = pos_start .+ first(element.LCS) * xinc'

    return ElementDisplacements(element, resolution, xinc, D, Dglobal, basepoints)
end

function ElementDisplacements(elements::AbstractVector{<:Asap.AbstractElement}, model::Asap.Model; resolution = 20)

    xstore = Vector{Float64}()
    ulocalstore = Vector{Matrix{Float64}}()
    uglobalstore = Vector{Matrix{Float64}}()
    basepointstore = Vector{Matrix{Float64}}()

    resolution = Int(round(resolution / length(elements)))

    element_loads = get_elemental_loads(model)

    for element in elements
        # Strip units from length (convert Quantity → Float64)
        L = ustrip(u"m", uconvert(u"m", element.length))

        xinc = collect(range(0, L, resolution))

        Dx, Dy, Dz = unodal(element; n = resolution)

        for load in element_loads[element.elementID]
            accumulatedisp!(load, xinc, Dy, Dz)
        end

        D = [Dx'; Dy'; Dz']
        Dglobal = hcat([sum(Δ .* element.LCS) for Δ in eachcol(D)]...)

        # Strip units from position (convert Quantity → Float64) for basepoints calculation
        pos_start = [ustrip(u"m", uconvert(u"m", p)) for p in element.nodeStart.position]
        basepoints = pos_start .+ first(element.LCS) * xinc'

        if isempty(xstore)
            xstore = [xstore; xinc]
        else
            xstore = [xstore; xstore[end] .+ xinc]
        end

        push!(ulocalstore, D)
        push!(uglobalstore, Dglobal)
        push!(basepointstore, basepoints)
    end

    return ElementDisplacements(elements[1], resolution, xstore, hcat(ulocalstore...), hcat(uglobalstore...), hcat(basepointstore...))
end

"""
    displacements(model::Asap.Model, increment)

Get the displacements of all elements in a model.

# Arguments
- `model::Model` - Structural model to analyze
- `increment` - Distance between sampling points (accepts Unitful length or Real in meters)
"""
function displacements(model::Asap.Model, increment)
    # Accept Unitful or Real - strip to meters internally
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)
    
    results = Vector{ElementDisplacements}()
    ids = groupbyid(model.elements)

    for id in ids
        elements = model.elements[id]
        L = sum(ustrip(u"m", el.length) for el in elements)
        n = max(Int(round(L / inc_m)), 2)

        push!(results, ElementDisplacements(elements, model; resolution = n))
    end

    return results
end