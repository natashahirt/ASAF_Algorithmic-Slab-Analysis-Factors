
struct GHnode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    reaction::Vector{Float64}
    u::Vector{Float64}
    displacement::Vector{Float64}
    id::String

    function GHnode(node::TrussNode)
        # Strip units from position (convert Quantity → Float64)
        position = [ustrip(u"m", uconvert(u"m", p)) for p in node.position]
        dof = [node.dof; [true, true, true]]
        nodeID = node.nodeID - 1
        # Strip units from reactions (N → Float64) and displacements (m → Float64)
        reaction = [ustrip(u"N", uconvert(u"N", r)) for r in node.reaction]
        u = [ustrip(u"m", uconvert(u"m", d)) for d in node.displacement]
        displacement = [ustrip(u"m", uconvert(u"m", d)) for d in node.displacement[1:3]]
        id = isnothing(node.id) ? "" : string(node.id)

        return new(position, dof, nodeID, reaction, u, displacement, id)
    end

    function GHnode(node::Node)
        # Strip units from position (convert Quantity → Float64)
        position = [ustrip(u"m", uconvert(u"m", p)) for p in node.position]
        dof = node.dof
        nodeID = node.nodeID - 1
        # Strip units from reactions and displacements
        # Reactions: DOFs 1-3 are forces (N), DOFs 4-6 are moments (N*m)
        # Displacements: DOFs 1-3 are translations (m), DOFs 4-6 are rotations (rad)
        reaction = Vector{Float64}(undef, length(node.reaction))
        for (i, r) in enumerate(node.reaction)
            if i <= 3
                reaction[i] = ustrip(u"N", uconvert(u"N", r))
            else
                reaction[i] = ustrip(u"N*m", uconvert(u"N*m", r))
            end
        end
        u = Vector{Float64}(undef, length(node.displacement))
        for (i, d) in enumerate(node.displacement)
            if i <= 3
                u[i] = ustrip(u"m", uconvert(u"m", d))
            else
                u[i] = ustrip(d)  # radians are dimensionless
            end
        end
        displacement = [ustrip(u"m", uconvert(u"m", d)) for d in node.displacement[1:3]]
        id = isnothing(node.id) ? "" : string(node.id)

        return new(position, dof, nodeID, reaction, u, displacement, id)
    end
end



struct GHsection
    E::Float64
    G::Float64
    A::Float64
    Ix::Float64
    Iy::Float64
    J::Float64

    function GHsection(section::TrussSection)
        # Strip units from section properties (convert Quantity → Float64)
        E = ustrip(u"Pa", uconvert(u"Pa", section.E))
        A = ustrip(u"m^2", uconvert(u"m^2", section.A))
        return new(E, 1.0, A, 1.0, 1.0, 1.0)
    end

    function GHsection(section::Section)
        # Strip units from section properties (convert Quantity → Float64)
        E = ustrip(u"Pa", uconvert(u"Pa", section.E))
        G = ustrip(u"Pa", uconvert(u"Pa", section.G))
        A = ustrip(u"m^2", uconvert(u"m^2", section.A))
        Ix = ustrip(u"m^4", uconvert(u"m^4", section.Ix))
        Iy = ustrip(u"m^4", uconvert(u"m^4", section.Iy))
        J = ustrip(u"m^4", uconvert(u"m^4", section.J))
        return new(E, G, A, Ix, Iy, J)
    end
end

const release2bool = Dict(
    Asap.FixedFixed => [false, false, false, false, false, false],
    Asap.FreeFixed => [true, true, true, false, false, false],
    Asap.FixedFree => [false, false, false, true, true, true],
    Asap.FreeFree => [true, true, true, true, true, true],
    Asap.Joist => [false, true, true, false, true, true]
)

function convert_release_to_bool(element::Element{T}) where {T}

    return release2bool[T]

end

struct GHelement
    iStart::Int64
    iEnd::Int64
    elementID::Int64
    section::GHsection
    release::Vector{Bool}
    psi::Float64
    localx::Vector{Float64}
    localy::Vector{Float64}
    localz::Vector{Float64}
    forces::Vector{Float64}
    axialforce::Float64
    id::String

    function GHelement(element::TrussElement)
        istart, iend = Asap.nodeids(element) .- 1
        elementID = element.elementID - 1
        section = GHsection(element.section)
        psi = element.Ψ
        lx, ly, lz = element.LCS
        id = isnothing(element.id) ? "" : string(element.id)
        forces = element.forces
        axialforce = forces[2]
        release = fill(true, 6)

        return new(istart, iend, elementID, section, release, psi, lx, ly, lz, forces, axialforce, id)
    end
    
    function GHelement(element::Element)
        istart, iend = Asap.nodeids(element) .- 1
        elementID = element.elementID - 1
        section = GHsection(element.section)
        psi = element.Ψ
        lx, ly, lz = element.LCS
        id = isnothing(element.id) ? "" : string(element.id)
        forces = element.forces
        axialforce = forces[7]
        release = convert_release_to_bool(element)

        return new(istart, iend, elementID, section, release, psi, lx, ly, lz, forces, axialforce, id)
    end
end

abstract type GHload end
function GHload end

struct GHnodeforce <: GHload
    value::Vector{Float64}
    id::String
    iNode::Int64
end

function GHload(load::NodeForce)
    # Strip units from load value (convert Quantity → Float64)
    i = load.node.nodeID - 1
    value = [ustrip(u"N", uconvert(u"N", v)) for v in load.value]
    id = isnothing(load.id) ? "" : string(load.id)

    return GHnodeforce(value, id, i)
end

struct GHnodemoment <: GHload
    value::Vector{Float64}
    id::String
    iNode::Int64
end

function GHload(load::NodeMoment)
    # Strip units from load value (convert Quantity → Float64)
    i = load.node.nodeID - 1
    value = [ustrip(u"N*m", uconvert(u"N*m", v)) for v in load.value]
    id = isnothing(load.id) ? "" : string(load.id)

    return GHnodemoment(value, id, i)
end

struct GHlineload <: GHload
    value::Vector{Float64}
    id::String
    iElement::Int64
end

function GHload(load::LineLoad)
    # Strip units from load value (convert Quantity → Float64)
    i = load.element.elementID - 1
    value = [ustrip(u"N/m", uconvert(u"N/m", v)) for v in load.value]
    id = isnothing(load.id) ? "" : string(load.id)

    return GHlineload(value, id, i)
end

struct GHpointload <: GHload
    value::Vector{Float64}
    id::String
    iElement::Int64
    x::Float64
end

function GHload(load::PointLoad)
    # Strip units from load value (convert Quantity → Float64)
    i = load.element.elementID - 1
    value = [ustrip(u"N", uconvert(u"N", v)) for v in load.value]
    id = isnothing(load.id) ? "" : string(load.id)
    x = load.position  # position is dimensionless Float64

    return GHpointload(value, id, i, x)
end

const loadtype2vectorindex = Dict(
    GHnodeforce => 1,
    GHnodemoment => 2,
    GHlineload => 3,
    GHpointload => 4
)

function categorize_loads(loads::Vector{<:GHload})
    nodeforces = Vector{GHnodeforce}()
    nodemoments = Vector{GHnodemoment}()
    lineloads = Vector{GHlineload}()
    pointloads = Vector{GHpointload}()

    load_collectors = [nodeforces, nodemoments, lineloads, pointloads]

    for load in loads
        i = loadtype2vectorindex[typeof(load)]
        push!(load_collectors[i], load)
    end

    return load_collectors
end

struct GHmodel
    nodes::Vector{GHnode}
    elements::Vector{GHelement}
    nodeforces::Vector{GHnodeforce}
    nodemoments::Vector{GHnodemoment}
    lineloads::Vector{GHlineload}
    pointloads::Vector{GHpointload}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}
    istart::Vector{Int64}
    iend::Vector{Int64}
    i_free_nodes::Vector{Int64}
    i_fixed_nodes::Vector{Int64}
end

"""
    GHmodel(model::Asap.AbstractModel)::GHmodel

Convert an Asap model into a condensed geometric format for interpoperability with other software.
"""
function GHmodel(model::Asap.AbstractModel)

    model.processed || (Asap.process!(model))

    nodes = GHnode.(model.nodes)
    elements = GHelement.(model.elements)
    loads = GHload.(model.loads)

    nodeforces, nodemoments, lineloads, pointloads = categorize_loads(loads)

    xyz = node_positions(model)

    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]

    istart = getproperty.(elements, :iStart)
    iend = getproperty.(elements, :iEnd)

    ifree = Vector{Int64}()
    ifixed = Vector{Int64}()

    dx = zero(x)
    dy = zero(y)
    dz = zero(z)

    for i in eachindex(nodes)

        if all(nodes[i].dof)
            push!(ifree, nodes[i].nodeID)
        else
            push!(ifixed, nodes[i].nodeID)
        end

        disp = nodes[i].displacement
        dx[i] = disp[1]
        dy[i] = disp[2]
        dz[i] = disp[3]
    end

    return GHmodel(
        nodes,
        elements,
        nodeforces,
        nodemoments,
        lineloads,
        pointloads,
        x,
        y,
        z,
        dx,
        dy,
        dz,
        istart,
        iend,
        ifree,
        ifixed
    )

end


"""
    GHsave(model::Asap.AbstractModel, filename::String)

Save an AsapModel as a .json file with GHmodel/GHnode/GHelement/GHload data structure. Add ".json" to the end of `filename` is optional.
"""
function GHsave(model::Asap.AbstractModel, filename::String)
    if filename[end-4:end] != ".json"
        filename *= ".json"
    end

    ghmodel = GHmodel(model)
    open(filename, "w") do f
        write(f, JSON.json(ghmodel))
    end
end

"""
    GHsave(model::GHmodel, filename::String)

Save a GHmodel as a .json file with GHmodel/GHnode/GHelement/GHload data structure. Add ".json" to the end of `filename` is optional.
"""
function GHsave(model::GHmodel, filename::String)
    if filename[end-4:end] != ".json"
        filename *= ".json"
    end

    open(filename, "w") do f
        write(f, JSON.json(model))
    end
end