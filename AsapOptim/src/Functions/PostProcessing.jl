"""
    Flocal(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
[2×1] vector of end element end forces in LCS
"""
function Flocal(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::AbstractOptParams)

    #Displacements w/r/t each element
    uEs = [u[id] for id in p.dofids]

    #End forces
    Rs .* Eks .* uEs
end

function ChainRulesCore.rrule(::typeof(Flocal), u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)

    F = Flocal(u, Eks, Rs, p)

    function Flocal_pullback(F̄)
        du = zero(u)
        dK = zero.(Eks)
        dR = zero.(Rs)
        ids = p.dofids

        for i in eachindex(ids)
            # F̄ ⋅ dF/dR
            ΔR = kron(Eks[i] * u[ids[i]], I(2)) * F̄[i]
            dR[i] = reshape(ΔR, 2, 6)

            #F̄ ⋅ dF/du
            du[ids[i]] += (Rs[i] * Eks[i])' * F̄[i]

            #F̄ ⋅ dF/dK
            ΔK = kron(u[ids[i]]', Rs[i])' * F̄[i]
            dK[i] = reshape(ΔK, 6, 6)
        end

        return (NoTangent(), du, dK, dR, NoTangent())

    end

    return F, Flocal_pullback
end

function ChainRulesCore.rrule(::typeof(Flocal), u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::FrameOptParams)

    F = Flocal(u, Eks, Rs, p)

    function Flocal_pullback(F̄)
        du = zero(u)
        dK = zero.(Eks)
        dR = zero.(Rs)
        ids = p.dofids

        for i in eachindex(ids)
            # F̄ ⋅ dF/dR
            ΔR = kron(Eks[i] * u[ids[i]], I(2)) * F̄[i]
            dR[i] = reshape(ΔR, 12, 12)

            #F̄ ⋅ dF/du
            du[ids[i]] += (Rs[i] * Eks[i])' * F̄[i]

            #F̄ ⋅ dF/dK
            ΔK = kron(u[ids[i]]', Rs[i])' * F̄[i]
            dK[i] = reshape(ΔK, 12, 12)
        end

        return (NoTangent(), du, dK, dR, NoTangent())

    end

    return F, Flocal_pullback
end

"""
    Flocal(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
[2×1] vector of end element end forces in LCS
"""
function axial_force(u::Vector{Float64}, Eks::Vector{Matrix{Float64}}, Rs::Vector{Matrix{Float64}}, p::TrussOptParams)
    fvecs = Flocal(u, Eks, Rs, p)
    getindex.(fvecs, 2)
end

"""
    σaxial(F::Vector{Float64}, A::Vector{Float64})

Get the absolute value of the axial stress experienced under force F and area A
"""
function σaxial(F::Vector{Float64}, A::Vector{Float64})
    stress = zero(F)

    for i in eachindex(F)
        if F[i] < 0
            stress[i] = -F[i] / A[i]
        else
            stress[i] = F[i] / A[i]
        end
    end
    

    return stress
end

"""
dg/dF = dg/dσ ⋅ dσ/dF = σ̄ ⋅ F / |F| / A
dg/dA = dg/dσ ⋅ dσ/dA = σ̄ ⋅ -F / A²
"""
function ChainRulesCore.rrule(::typeof(σaxial), F::Vector{Float64}, A::Vector{Float64})
    stress = σaxial(F, A)

    function σaxial_pullback(σ̄)
        s = sign.(F)
        return (NoTangent(), σ̄ .* s ./ A, -σ̄ .* F ./ A.^2)
    end

    return stress, σaxial_pullback
end

"""
    axialforce(t::TrussResults, p::TrussOptParams)
Axial forces in truss structure
"""
function axial_force(t::TrussResults, p::TrussOptParams)
    axial_force(t.U, t.K, t.R, p)
end

"""
    axial_stress(t::TrussResults, p::TrussOptParams)
Axial stresses in truss structure
"""
function axial_stress(t::TrussResults, p::TrussOptParams)
    F = axial_force(t.U, t.K, t.R, p)
    σaxial(F, t.A)
end

"""
    axialforce(t::NetworkResults, p::NetworkOptParams)
Axial forces in FDM network
"""
function axial_force(t::NetworkResults, p::NetworkOptParams)
    norm.(eachrow(p.C * [t.X t.Y t.Z])) .* t.Q
end

"""
    updatemodel(p::TrussOptParams, u::Vector{Float64})

Generate a new structural model from the results of an optimization
"""
function updatemodel(p::TrussOptParams, u::Vector{Float64})
    
    #final values
    X = add_values(p.X, p.indexer.iX, u[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, u[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, u[p.indexer.iZg] .* p.indexer.fZ)
    A = replace_values(p.A, p.indexer.iA, u[p.indexer.iAg] .* p.indexer.fA)

    #new model
    nodes = Vector{TrussNode}()
    elements = Vector{TrussElement}()
    loads = Vector{NodeForce}()

    #new nodes
    for (node, x, y, z) in zip(p.model.nodes, X, Y, Z)
        newnode = TrussNode([x, y, z], node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for (id, e, a, el) in zip(p.nodeids, p.E, A, p.model.elements)
        newelement = TrussElement(nodes[id]..., TrussSection(a, e))
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.model.loads
        newload = NodeForce(nodes[load.node.nodeID], load.value)
        newload.id = load.id
        push!(loads, newload)
    end

    model = TrussModel(nodes, elements, loads)
    Asap.solve!(model)

    return model

end

"""
    updatenetwork(p::TrussOptParams, u::Vector{Float64})

Generate a new FDM network from the results of an optimization
"""
function updatenetwork(p::NetworkOptParams, u::Vector{Float64})
    
    #final values
    X = add_values(p.X, p.indexer.iX, u[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, u[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, u[p.indexer.iZg] .* p.indexer.fZ)
    Q = replace_values(p.q, p.indexer.iQ, u[p.indexer.iQg] .* p.indexer.fQ)

    #new model
    nodes = Vector{FDMnode}()
    elements = Vector{FDMelement}()
    loads = Vector{FDMload}()

    #new nodes
    for (node, x, y, z) in zip(p.network.nodes, X, Y, Z)
        newnode = FDMnode([x, y, z], node.dof)
        newnode.id = node.id
        push!(nodes, newnode)
    end

    #new elements
    for (q, el) in zip(Q, p.network.elements)
        newelement = FDMelement(nodes, el.iStart, el.iEnd, q)
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.network.loads
        newload = FDMload(nodes[load.point.nodeID], load.force)
        push!(loads, newload)
    end

    network = Network(nodes, elements, loads)
    Asap.solve!(network)

    return network

end

"""
    updatemodel(p::TrussOptParams, u::Vector{Float64})

Generate a new structural model from the results of an optimization
"""
function updatemodel(p::FrameOptParams, u::Vector{Float64})
    
    #final values
    X = add_values(p.X, p.indexer.iX, u[p.indexer.iXg] .* p.indexer.fX)
    Y = add_values(p.Y, p.indexer.iY, u[p.indexer.iYg] .* p.indexer.fY)
    Z = add_values(p.Z, p.indexer.iZ, u[p.indexer.iZg] .* p.indexer.fZ)

    A = replace_values(p.A, p.indexer.iA, u[p.indexer.iAg] .* p.indexer.fA)
    Ix = replace_values(p.Ix, p.indexer.iIx, u[p.indexer.iIxg] .* p.indexer.fIx)
    Iy = replace_values(p.Iy, p.indexer.iIy, u[p.indexer.iIyg] .* p.indexer.fIy)
    J = replace_values(p.J, p.indexer.iJ, u[p.indexer.iJg] .* p.indexer.fJ)

    #new model
    nodes = Vector{Node}()
    elements = Vector{Element}()
    loads = Vector{Asap.AbstractLoad}()

    #new nodes
    for (node, x, y, z) in zip(p.model.nodes, X, Y, Z)
        newnode = Node([x, y, z], node.dof, node.id)
        push!(nodes, newnode)
    end

    #new elements
    for (id, el, a, ix, iy, j) in zip(p.nodeids, p.model.elements, A, Ix, Iy, J)
        section = deepcopy(el.section)

        newsection = Section(a, section.E, section.G, ix, iy, j)
        newelement = Element(nodes[id]..., newsection)
        newelement.id = el.id
        push!(elements, newelement)
    end

    #new loads
    for load in p.model.loads

        if typeof(load) <: Asap.NodeLoad
            newload = typeof(load)(nodes[load.node.nodeID], load.value)
            newload.id = load.id
            push!(loads, newload)
        else
            if typeof(load) == Asap.PointLoad
                newload = PointLoad(elements[load.element.elementID], load.position, load.value)
                newload.id = load.id
                push!(loads, newload)
            else
                newload = typeof(load)(elements[load.element.elementID], load.value)
                newload.id = load.id
                push!(loads, newload)
            end
        end
    end

    model = Model(nodes, elements, loads)
    Asap.solve!(model)

    return model

end

