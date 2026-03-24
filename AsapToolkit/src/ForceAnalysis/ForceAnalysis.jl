"""
Accumlate the internal forces cause by a given load to the current element
"""
function accumulate_force!(load::LineLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", load.element.length))

    # distributed load magnitudes in LCS (strip units from load.value)
    value_stripped = [ustrip(u"N/m", uconvert(u"N/m", v)) for v in load.value]
    wx, wy, wz = (R * value_stripped) .* [1, -1, -1]

    P .+= [PLine(wx, L, x) for x in xvals]

    My .+= [MLine(load.element, wy, L, x) for x in xvals]
    Vy .+= [VLine(load.element, wy, L, x) for x in xvals]

    Mz .+= [MLine(load.element, wz, L, x) for x in xvals]
    Vz .+= [VLine(load.element, wz, L, x) for x in xvals]
end

"""
Accumlate the internal forces cause by a given load to the current element
"""
function accumulate_force!(load::PointLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", load.element.length))
    frac = load.position

    # distributed load magnitudes in LCS (strip units from load.value)
    value_stripped = [ustrip(u"N", uconvert(u"N", v)) for v in load.value]
    px, py, pz = (R * value_stripped) .* [1, -1, -1]

    P .+= [PPoint(px, L, x, frac) for x in xvals]

    My .+= [MPoint(load.element, py, L, x, frac) for x in xvals]
    Vy .+= [VPoint(load.element, py, L, x, frac) for x in xvals]

    Mz .+= [MPoint(load.element, pz, L, x, frac) for x in xvals]
    Vz .+= [VPoint(load.element, pz, L, x, frac) for x in xvals]
end

"""
Accumulate the internal forces caused by a GravityLoad (self-weight) to the current element.
GravityLoad applies a distributed load equal to ρ * A * g in the global -Z direction.
"""
function accumulate_force!(load::Asap.GravityLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    element = load.element
    R = element.R[1:3, 1:3]
    # Strip units from length
    L = ustrip(u"m", uconvert(u"m", element.length))

    # Self-weight per unit length: w = ρ * A * g (in N/m)
    ρ = ustrip(u"kg/m^3", element.section.ρ)
    A = ustrip(u"m^2", element.section.A)
    g = ustrip(u"m/s^2", load.factor)
    w_mag = ρ * A * g  # N/m

    # Global load vector: self-weight acts in -Z direction
    w_global = [0.0, 0.0, -w_mag]

    # Transform to local coordinate system
    wx, wy, wz = (R * w_global) .* [1, -1, -1]

    P .+= [PLine(wx, L, x) for x in xvals]

    My .+= [MLine(element, wy, L, x) for x in xvals]
    Vy .+= [VLine(element, wy, L, x) for x in xvals]

    Mz .+= [MLine(element, wz, L, x) for x in xvals]
    Vz .+= [VLine(element, wz, L, x) for x in xvals]
end

"""
Store the internal force results for a given element
"""
struct InternalForces
    element::Element
    resolution::Integer
    x::Vector{Float64}
    P::Vector{Float64}
    My::Vector{Float64}
    Vy::Vector{Float64}
    Mz::Vector{Float64}
    Vz::Vector{Float64}
end


"""
Internal force sampling for an element 
"""
function InternalForces(element::Asap.AbstractElement, model::Asap.Model; resolution = 20)
    
    #beam information
    dofs = etype2DOF[typeof(element)]
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", element.length))

    #discretization
    xinc = collect(range(0, L, resolution))

    #end node information
    # Strip units from displacements for matrix multiplication (R and K are Float64)
    uglobal = [Asap.to_displacement_vec(element.nodeStart.displacement); Asap.to_displacement_vec(element.nodeEnd.displacement)]
    
    # end forces that are relevant to the given element/release condition
    Flocal = (element.R * element.K * uglobal) .* dofs

    # shear/moment acting at the *starting* point of an element in LCS
    Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

    # initialize internal force vectors
    P = repeat([Pstart], resolution)
    My = Vystart .* xinc .- Mystart
    Vy = zero(My) .+ Vystart
    Mz = Vzstart .* xinc .- Mzstart
    Vz = zero(Mz) .+ Vzstart

    element_loads = get_elemental_loads(model)

    # accumulate loads
    for load in element_loads[element.elementID]
        accumulate_force!(load,
            xinc,
            P,
            My,
            Vy,
            Mz,
            Vz)  
    end

    return InternalForces(element, resolution, xinc, P, My, Vy, Mz, Vz)
end

"""
    InternalForces(element::Element, loads::Vector{<:ElementLoad}; resolution = 20)

Get internal force results for a given element from a set of loads
"""
function InternalForces(element::Asap.AbstractElement, loads::AbstractVector{<:Asap.AbstractLoad}; resolution = 20)
    
    #beam information
    dofs = etype2DOF[typeof(element)]
    # Strip units from length (convert Quantity → Float64)
    L = ustrip(u"m", uconvert(u"m", element.length))

    #discretization
    xinc = collect(range(0, L, resolution))

    #end node information
    # Strip units from displacements for matrix multiplication (R and K are Float64)
    uglobal = [Asap.to_displacement_vec(element.nodeStart.displacement); Asap.to_displacement_vec(element.nodeEnd.displacement)]
    
    # end forces that are relevant to the given element/release condition
    Flocal = (element.R * element.K * uglobal) .* dofs

    # shear/moment acting at the *starting* point of an element in LCS
    Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

    # initialize internal force vectors
    P = repeat([Pstart], resolution)
    My = Vystart .* xinc .- Mystart
    Vy = zero(My) .+ Vystart
    Mz = Vzstart .* xinc .- Mzstart
    Vz = zero(Mz) .+ Vzstart

    # accumulate loads
    for load in loads
        accumulate_force!(load,
            xinc,
            P,
            My,
            Vy,
            Mz,
            Vz)  
    end

    return InternalForces(element, resolution, xinc, P, My, Vy, Mz, Vz)
end


function get_elemental_loads(model::Model)

    element_to_loads = [AbstractLoad[] for _ in 1:model.nElements]
    for load in model.loads
        hasproperty(load, :element) || continue
        push!(element_to_loads[load.element.elementID], load)
    end

    return element_to_loads

end

"""
    InternalForces(element::Vector{<:FrameElement}, model::Model; resolution = 20)

Get internal force results for a group of ordered elements that form a single physical element
"""
function InternalForces(elements::AbstractVector{<:Asap.AbstractElement}, model::Asap.Model; resolution = 20)
    
    xstore = Vector{Float64}()
    pstore = Vector{Float64}()
    mystore = Vector{Float64}()
    vystore = Vector{Float64}()
    mzstore = Vector{Float64}()
    vzstore = Vector{Float64}()

    resolution = Int(round(resolution / length(elements)))

    element_loads_map = get_elemental_loads(model)

    #beam information
    for element in elements
        loadids = element_loads_map[element.elementID]

        dofs = etype2DOF[typeof(element)]
        # Strip units from length (convert Quantity → Float64)
        L = ustrip(u"m", uconvert(u"m", element.length))

        #discretization
        xinc = collect(range(0, L, max(resolution, 2)))

        #end node information
        uglobal = [element.nodeStart.displacement; element.nodeEnd.displacement]
        
        # end forces that are relevant to the given element/release condition
        Flocal = (element.R * element.K * uglobal) .* dofs

        # shear/moment acting at the *starting* point of an element in LCS
        Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

        # initialize internal force vectors
        P = repeat([Pstart], resolution)
        My = Vystart .* xinc .- Mystart
        Vy = zero(My) .+ Vystart
        Mz = Vzstart .* xinc .- Mzstart
        Vz = zero(Mz) .+ Vzstart

        # accumulate loads
        for load in loadids
            accumulate_force!(load,
                xinc,
                P,
                My,
                Vy,
                Mz,
                Vz)  
        end

        if isempty(xstore)
            xstore = [xstore; xinc]
        else
            xstore = [xstore; xstore[end] .+ xinc]
        end

        pstore = [pstore; P]
        mystore = [mystore; My]
        vystore = [vystore; Vy]
        mzstore = [mzstore; Mz]
        vzstore = [vzstore; Vz]

    end

    return InternalForces(elements[1], resolution, xstore, pstore, mystore, vystore, mzstore, vzstore)
end



"""
    forces(model::Asap.Model, increment::Real)

Get the internal forces of all elements in a model.

# Arguments
- `model::Model` - Structural model whose elements are analyzed
- `increment` - Distance between sampling points (accepts Unitful length or Real in meters)
"""
function InternalForces(model::Asap.Model, increment)
    # Accept Unitful or Real - strip to meters internally
    inc_m = increment isa Unitful.Quantity ? ustrip(u"m", increment) : Float64(increment)

    results = Vector{InternalForces}()
    element_loads_map = get_elemental_loads(model)

    for element in model.elements
        L = ustrip(u"m", element.length)
        n = max(Int(round(L / inc_m)), 2)

        push!(results, InternalForces(element, element_loads_map[element.elementID]; resolution = n))
    end

    return results
end

struct ForceEnvelopes
    element::Element
    resolution::Integer
    x::Vector{Float64}
    Plow::Vector{Float64}
    Phigh::Vector{Float64}
    Mylow::Vector{Float64}
    Myhigh::Vector{Float64}
    Vylow::Vector{Float64}
    Vyhigh::Vector{Float64}
    Mzlow::Vector{Float64}
    Mzhigh::Vector{Float64}
    Vzlow::Vector{Float64}
    Vzhigh::Vector{Float64}
end

"""
    load_envelopes(model::Model, loads::Vector{Vector{<:AbstractLoad}})

Get the high/low internal forces for a series of external loads
"""
function load_envelopes(model::Model, loads::Vector{Vector{<:Asap.AbstractLoad}}, increment::Real)
    #
    envelopes = Vector{ForceEnvelopes}()

    #collector of force results
    forceresults = Vector{Vector{InternalForces}}()

    # perform analysis
    for load in loads
        solve!(model, load)
        push!(forceresults, forces(model, increment))
    end

    # number of actual elements
    n = length(first(forceresults))

    # for each element
    for i = 1:n
        e = first(forceresults)[i].element
        res = first(forceresults)[i].resolution
        x = first(forceresults)[i].x

        P = hcat(getproperty.(getindex.(forceresults, i), :P)...)
        My = hcat(getproperty.(getindex.(forceresults, i), :My)...)
        Vy = hcat(getproperty.(getindex.(forceresults, i), :Vy)...)
        Mz = hcat(getproperty.(getindex.(forceresults, i), :Mz)...)
        Vz = hcat(getproperty.(getindex.(forceresults, i), :Vz)...)

        Prange = extrema.(eachrow(P))
        Myrange = extrema.(eachrow(My))
        Vyrange = extrema.(eachrow(Vy))
        Mzrange = extrema.(eachrow(Mz))
        Vzrange = extrema.(eachrow(Vz))

        envelope = ForceEnvelopes(e,
            res,
            x,
            getindex.(Prange, 1),
            getindex.(Prange, 2),
            getindex.(Myrange, 1),
            getindex.(Myrange, 2),
            getindex.(Vyrange, 1),
            getindex.(Vyrange, 2),
            getindex.(Mzrange, 1),
            getindex.(Mzrange, 2),
            getindex.(Vzrange, 1),
            getindex.(Vzrange, 2))

        push!(envelopes, envelope)
    end

    return envelopes
end