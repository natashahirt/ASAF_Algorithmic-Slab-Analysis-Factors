struct ModelGeo <: AbstractGeo
    nodes::Vector{Vector{Float64}}
    nodes_xy::Vector{Vector{Float64}}
    disp::Vector{Vector{Float64}}
    disp_xy::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    indices_flat::Vector{Int64}
    P::Vector{Float64}
    max_abs_P::Float64
    Vy::Vector{Float64}
    max_abs_Vy::Float64
    Vz::Vector{Float64}
    max_abs_Vz::Float64
    Tx::Vector{Float64}
    max_abs_Tx::Float64
    My::Vector{Float64}
    max_abs_My::Float64
    Mz::Vector{Float64}
    max_abs_Mz::Float64
    areas::Vector{Float64}
    max_area::Float64
    Ix::Vector{Float64}
    max_Ix::Float64
    Iy::Vector{Float64}
    max_Iy::Float64
    J::Vector{Float64}
    max_J::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    element_vectors_xy::Vector{Vector{Float64}}

    function ModelGeo(model::Model)
        # Strip units from positions (convert Quantity → Float64)
        nodes = [[ustrip(u"m", uconvert(u"m", p)) for p in node.position] for node in model.nodes]
        nodes_xy = [node[1:2] for node in nodes]

        # Strip units from displacements: first 3 DOFs are translations (m)
        disp = [Asap.to_displacement_vec(node.displacement)[1:3] for node in model.nodes]
        disp_xy = [d[1:2] for d in disp]

        indices = Asap.nodeids.(model.elements)
        indices_flat = vcat(indices...)

        element_forces = getproperty.(model.elements, :forces)

        P = vcat([forces[[1,7]] .* [-1, 1] for forces in element_forces]...)
        max_abs_P = maximum(abs.(P))

        Vy = vcat([forces[[2,8]] .* [-1, 1] for forces in element_forces]...)
        max_abs_Vy = maximum(abs.(Vy))

        Vz = vcat([forces[[3,9]] .* [-1, 1] for forces in element_forces]...)
        max_abs_Vz = maximum(abs.(Vz))

        Tx = vcat([forces[[4,10]] .* [-1, 1] for forces in element_forces]...)
        max_abs_Tx = maximum(abs.(Tx))

        My = vcat([forces[[5,11]] .* [-1, 1] for forces in element_forces]...)
        max_abs_My = maximum(abs.(My))

        Mz = vcat([forces[[6,12]] .* [-1, 1] for forces in element_forces]...)
        max_abs_Mz = maximum(abs.(Mz))

        sections = getproperty.(model.elements, :section)

        # Strip units from section properties (convert Quantity → Float64)
        areas = [ustrip(u"m^2", uconvert(u"m^2", s.A)) for s in sections]
        max_area = maximum(areas)

        Ix = [ustrip(u"m^4", uconvert(u"m^4", s.Ix)) for s in sections]
        max_Ix = maximum(Ix)

        Iy = [ustrip(u"m^4", uconvert(u"m^4", s.Iy)) for s in sections]
        max_Iy = maximum(Iy)

        J = [ustrip(u"m^4", uconvert(u"m^4", s.J)) for s in sections]
        max_J = maximum(J)

        element_vectors = Asap.local_x.(model.elements)
        element_vectors_xy = [evec[1:2] for evec in element_vectors]

        # Strip units from lengths (convert Quantity → Float64)
        lengths = [ustrip(u"m", uconvert(u"m", e.length)) for e in model.elements]

        return new(
            nodes,
            nodes_xy,
            disp,
            disp_xy,
            indices,
            indices_flat,
            P,
            max_abs_P,
            Vy,
            max_abs_Vy,
            Vz,
            max_abs_Vz,
            Tx,
            max_abs_Tx,
            My,
            max_abs_My,
            Mz,
            max_abs_Mz,
            areas,
            max_area,
            Ix,
            max_Ix,
            Iy,
            max_Iy,
            J,
            max_J,
            lengths,
            element_vectors,
            element_vectors_xy
        )
    end

end