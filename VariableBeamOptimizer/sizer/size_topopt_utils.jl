struct CellData
    area::Float64
    distances::Vector{Float64}
    beam_idxs::Vector{Int}
    beam_t::Vector{Float64}
end

function suppress_if_small(x, variable; sharpness::Real=20.0)
    threshold_dict = Dict(:h => 0.01, :w => 0.01, :tw => 0.001, :tf => 0.001)
    return x * (0.5 * (1 + tanh(sharpness * (x / threshold_dict[variable] - 1)))) + 1e-10  # ∈ [1e-10, x]
end    

function get_volume(minimizers::Union{Vector{Vector{Float64}}, Vector{Vector}}, model::Asap.Model{Element, Asap.AbstractLoad})
    beam_elements = model.elements[:beam]
    return sum([I_symm(minimizers[i]...).A * (beam_elements[i].length * 39.3701) for i in 1:lastindex(minimizers)])
end

function get_volume(minimizers::Vector{Vector{Float64}}, elements::Vector{Asap.Element})
    return sum([I_symm(minimizers[i]...).A * (elements[i].length * 39.3701) for i in 1:lastindex(minimizers)])
end

function get_volume(minimizers::Vector{Float64}, elements::Vector{Asap.Element})
    return sum([minimizers[i] * (elements[i].length * 39.3701) for i in 1:lastindex(minimizers)])
end

function get_norm_mass(volume, area)
    return volume * convert_to_m[:in]^3 * ρ_STEEL / area
end

function compute_beam_stiffnesses(h, w, tw, tf, l)
    # Compute the moment of inertia for each beam and return the stiffness
    return [Ix_I_symm(h[i], w[i], tw[i], tf[i]) / l[i]^3 for i in 1:length(h)]
end

function redistribute_cell(cell::CellData, stiffness_beams, factored_w_slab, factored_w_applied, slab_depth)
    d = cell.distances
    i = cell.beam_idxs
    t = cell.beam_t
    
    s = [1 / (d[j]^3 + 1e-6) + stiffness_beams[i[j]] for j in 1:length(i)]
    
    w = s ./ sum(s) # weights for each beam
    w_load = w * cell.area * (factored_w_slab * slab_depth + factored_w_applied) # loads per beam
    return i, w_load, t
end

function collect_load_triplets(cell_data, stiffness_beams, factored_w_slab, factored_w_applied, slab_depth)
    triplets = Tuple{Int, Float64, Float64}[]
    for cell in cell_data
        i_vec, w_vec, t_vec = redistribute_cell(cell[2], stiffness_beams, factored_w_slab, factored_w_applied, slab_depth)
        append!(triplets, zip(i_vec, t_vec, w_vec))
    end
    return triplets
end

function sum_vector(vectors::Vector{<:AbstractVector{<:Real}})
    reduce((a, b) -> a .+ b, vectors)
end

function compute_M_V_from_triplets(triplets, model; resolution=200)
    elements = model.elements[:beam]
    n_beams = length(elements)

    # === Precompute scalar/vector arrays for each beam ===
    Ls = [el.length for el in elements]
    Rs = [el.R[1:3, 1:3] for el in elements]
    disps = [vcat(el.nodeStart.displacement, el.nodeEnd.displacement) for el in elements]
    r2dofs = [AsapToolkit.release2DOF[AsapToolkit.get_release(el)] for el in elements]
    Flocals = [(el.R * el.K * disp) .* r2dof for (el, disp, r2dof) in zip(elements, disps, r2dofs)]
    Vystarts = [F[2] for F in Flocals]
    Mystarts = [F[6] for F in Flocals]
    xrel = collect(0:(resolution - 1)) ./ (resolution - 1)  # relative positions
    xincs = [xrel .* Ls[i] for i in 1:n_beams]  # scale by each beam's length
    releases = [AsapToolkit.get_release(el) for el in elements]
    moment_fns = [AsapToolkit.MPointLoad[r] for r in releases]
    shear_fns = [AsapToolkit.VPointLoad[r] for r in releases]

    My_all = [
        Vystarts[i] .* xincs[i] .- Mystarts[i] .+ sum_vector([moment_fns[i].((Rs[i] * [0.0, 0.0, -w])[2], Ls[i], xincs[i], t) for (j, t, w) in triplets if j == i])
        for i in 1:n_beams
    ]
    
    Vy_all = [
        fill(Vystarts[i], resolution) .+ sum_vector([shear_fns[i].((Rs[i] * [0.0, 0.0, -w])[2], Ls[i], xincs[i], t) for (j, t, w) in triplets if j == i])
        for i in 1:n_beams
    ]

    max_Ms = [maximum(abs.(My)) for My in My_all]
    max_Vs = [maximum(abs.(Vy)) for Vy in Vy_all]

    return vcat(max_Ms, max_Vs)
end

Zygote.@nograd function beam_demand_from_stiffness(stiffness_beams, cell_data, model, factored_w_slab, factored_w_applied, slab_depth)
    triplets = collect_load_triplets(cell_data, stiffness_beams, factored_w_slab, factored_w_applied, slab_depth)
    return compute_M_V_from_triplets(triplets, model)
end

Zygote.@nograd function solve_indeterminate_system(h, w, tw, tf, model, cell_data, factored_w_slab, factored_w_applied, slab_depth;
                                                   perimeter_indices::Vector{Int}=Int[], facade_line_load::Float64=0.0)
    
    model.loads = Asap.AbstractLoad[]
    beam_elements = model.elements[:beam]
    l = [element.length for element in beam_elements]

    stiffness_beams = compute_beam_stiffnesses(h, w, tw, tf, l)
    triplets = collect_load_triplets(cell_data, stiffness_beams, factored_w_slab, factored_w_applied, slab_depth)
    load_dict = Dict{Int, Vector{Asap.AbstractLoad}}()

    for (i, t, w) in triplets
        load = Asap.PointLoad(beam_elements[i], t, [0,0,-w])
        push!(model.loads, load)
        if !haskey(load_dict, i)
            load_dict[i] = Asap.AbstractLoad[]
        end
        push!(load_dict[i], load)
    end

    # Add perimeter facade line loads when requested.
    if !isempty(perimeter_indices) && !iszero(facade_line_load)
        for i in perimeter_indices
            line_load = LineLoad(beam_elements[i], [0,0,-facade_line_load])
            push!(model.loads, line_load)
            if !haskey(load_dict, i)
                load_dict[i] = Asap.AbstractLoad[]
            end
            push!(load_dict[i], line_load)
        end
    end

    for (i, element) in enumerate(beam_elements)
        I_symm_section = I_symm(h[i], w[i], tw[i], tf[i])
        element.section = Section(I_symm_section.A, steel_ksi.E, steel_ksi.G, I_symm_section.Ix, I_symm_section.Iy, I_symm_section.J)
    end

    Asap.solve!(model, reprocess=true)

    max_Ms = zeros(length(beam_elements))
    max_Vs = zeros(length(beam_elements))
    max_δs = zeros(length(beam_elements))
    
    for (i, el) in enumerate(beam_elements)
        internal_forces = InternalForces(el, get(load_dict, i, Asap.AbstractLoad[]); resolution=200)
        max_Ms[i] = maximum(abs.(internal_forces.My))
        max_Vs[i] = maximum(abs.(internal_forces.Vy))
        disp = ElementDisplacements(el, get(load_dict, i, Asap.AbstractLoad[]), resolution=200)
        max_δs[i] = maximum(abs.(disp.ulocal[2, :]))
    end

    return vcat(max_Ms, max_Vs, max_δs)

end