"""
    get_i_check(n_samples::Int64, n_sections::Int64)

takes the number of samples you've taken and the number of sections you want to analyze, and
returns n_sections number of evenly distributed indices in that range.

# Arguments
- n_sections::Int64 : number of sections that are desired
- beam_forces::InternalForces : need this to heuristically determine optimal section positions
"""
function get_i_check(n_sections::Int64, beam_forces::InternalForces; distribution::Symbol=:heuristic)

    @assert distribution in (:heuristic, :even)

    if n_sections < 2 # defer to max and min

        return Int64[]

    end
    
    n_samples = length(beam_forces.x)
    i_check = [1, n_samples]

    if distribution == :even

        n_sections -= 1
        increment = (n_samples-(n_samples%n_sections))/n_sections # given the number of sections we want to optimize, what's our increment
        
        append!(i_check, [Int(floor(i*increment)) for i in 1:n_sections][1:end-1])
    
    elseif distribution == :heuristic

        n_interior = n_sections - 2
        
        extreme_i_M = get_extreme_i(beam_forces.My)
        extreme_i_V = get_extreme_i(beam_forces.Vy)
        zero_i_M = get_roots_i(beam_forces.My)
        zero_i_V = get_roots_i(beam_forces.Vy)

        critical_i = collect(Set([extreme_i_M; extreme_i_V; zero_i_M; zero_i_V]))

        i_interior = [critical_i[i] for i in 1:min(n_interior, length(critical_i))]

        append!(i_check, i_interior)

        if length(i_check) < n_sections

            while length(i_check) < n_sections

                bin_ranges = get_bins(i_check)
                biggest_bin = bin_ranges[1]
                new_i = Int(floor((biggest_bin[1]+biggest_bin[2])/2))
                push!(i_check, new_i)

            end

        end

    end

    i_check = sort(i_check)

    return i_check

end

"""
    get_bins(i_check::Vector{Int})

given a list of indices (x1, x2, x3...xn) get the list of tuples that define the segments
(or bins) in between them (e.g. [(x1,x2),(x2,x3),...,(x(n-1),xn)])

# Arguments
- i_check::Vector{Int} : list of indices that we're getting the bins for
"""
function get_bins(i_check::Vector{Int})

    i_check = sort(i_check)
    bin_ranges = [(i_check[i], i_check[i+1]) for i in 1:lastindex(i_check)-1]
    bin_lengths = [bin[2]-bin[1] for bin in bin_ranges]
    # sort by biggest range first
    p = reverse(sortperm(bin_lengths))
    bin_ranges = bin_ranges[p]

    return bin_ranges

end

"""
    get_extreme_i(vector::Vector{T}) where T <: Real

iteratively find the extrema (i.e. where the second derivative changes sign) for a list of points

# Arguments
- vector::Vector{T} : vector of values we're trying to find the extrema for
"""
function get_extreme_i(vector::Vector{T}) where T <: Real

    extreme_i = []
    sign_i = sign(vector[1] - vector[2])

    for i in 2:lastindex(vector)-1

        new_sign = sign(vector[i] - vector[i+1])

        if sign_i * new_sign < 0

            push!(extreme_i, i)
            sign_i = new_sign

        end

    end

    return extreme_i

end

"""
    get_roots_i(vector::Vector{T}) where T <: Real

iteratively find the roots (i.e. where the series intersects the x axis) for a list of points

# Arguments
- vector::Vector{T} : vector of values we're trying to find the extrema for
"""
function get_roots_i(vector::Vector{T}) where T <: Real

    zero_i = []
    sign_i = sign(vector[1])

    for i in 2:lastindex(vector)-1

        new_sign = sign(vector[i])

        if sign_i * new_sign < 0

            push!(zero_i, i)
            sign_i = new_sign

        end

    end

    return zero_i

end


"""
    get_closest_i(interp_points::Vector{T}, target_x::T)

get the index of of the closest item in a list to the target value. Rounds down.

# Arguments
- interp_points::Vector{T} : list of values we're comparing the target to
- target_x::T : value that we're trying to get the closest index for
"""
function get_closest_i(interp_points::Vector{T}, target_x::T) where T <: Real

    if target_x > interp_points[end]

        return lastindex(interp_points)

    elseif target_x < interp_points[1]

        return 1
    
    else

        low = 1
        high = lastindex(interp_points)
        closest_i = 1
        closest_value = interp_points[closest_i]

        while low <= high
            
            mid_i = Int(floor((high+low)/2))
            midpoint = interp_points[mid_i]
    
            # Update closest if this value is closer to the target
            if abs(midpoint - target_x) < abs(closest_value - target_x)

                closest_value = midpoint
                closest_i = mid_i

            end
    
            # Move either low or high pointer
            if midpoint < target_x

                low = mid_i + 1

            else
                
                high = mid_i - 1

            end
    
        end

        return closest_i

    end

end


"""
    get_element_deflection(params::FrameOptParams, section::I_symm)

this is the CORRECT AND DIFFERENTIABLE deflection calculator I REST MY CASE
"""
function get_element_deflection(params::FrameOptParams, vars::Vector{Float64}; material::Union{Nothing,AbstractMaterial,Asap.Material}=steel_ksi)
    # Convert section properties from imperial (in) to metric (cm) units

    # Define variables to hold values
    A = A_I_symm(vars...) # in² 
    Ix = Ix_I_symm(vars...) # in⁴
    Iy = Iy_I_symm(vars...) # in⁴
    J = J_I_symm(vars...) # in⁴

    param_vars = collect(params.values)

    for i in 1:Int64(length(params.values)/4)
        param_vars = replace_values(param_vars, [4i-3,4i-2,4i-1,4i], [A, Ix, Iy, J])
    end

    results = AsapOptim.solve_frame_Pf(param_vars, params, dead_load = material.ρ) # yes dead load
    #results = AsapOptim.solve_frame(param_vars, params) # no dead load

    δ_local = results.U[3:6:end] # vertical deflection, in

    return δ_local # in

end

# Define a function to check if a load is of type PointLoad with a specific subtype of Asap.Release
function is_pointload(load)
    return load isa PointLoad{R} where R <: Asap.Release
end

function is_lineload(load)
    return load isa LineLoad{R} where R <: Asap.Release
end

function get_frameoptparams(sizing_params::SlabSizingParams, beam_element::Element, element_loads::Vector{Asap.AbstractLoad})

    loads = Asap.AbstractLoad[]
    length = beam_element.length
    all_load_values = [
        hasproperty(load, :loadID) ? get_new_load_value(sizing_params, getproperty(load, :loadID), factored=false) : load.value
        for load in element_loads
    ]

    section = beam_element.section
    E = section.E
    Ix = section.Ix

    if !isempty(all_load_values)
        w_estimated = sum(all_load_values) / (length)  # Assumixng all_load_values is the total load
    else
        w_estimated = 0.
    end

    if isempty(all_load_values)

        n_nodes = 100

        positions = collect(0:length/(n_nodes + 1):length) # in

        nodes = [Node([pos,0.,0.],Bool[1,1,1,0,0,0]) for pos in positions] # singular exception can be removed by replacing with Bool[1,0,0,0,1,0]    
        nodes[1].dof = Bool[0,0,0,1,1,1]
        nodes[end].dof = Bool[1,0,0,1,1,1]
        elements = [Element(nodes[n-1],nodes[n], section) for n in 2:lastindex(nodes)]

    elseif is_pointload(element_loads[1])

        all_positions = [load.position .* length for load in element_loads] # in
        
        positions = Float64[]
        load_values = Vector{Vector{Float64}}()

        # make sure duplicate positions aren't counted twice
        for i in 1:lastindex(all_positions)

            if all_positions[i] ∉ positions

                push!(positions, all_positions[i])
                push!(load_values, all_load_values[i])

            else

                indices = findall(x -> x == all_positions[i], positions)

                for index in indices

                    load_values[index] += all_load_values[i]

                end

            end

        end

        load_values = [[[0.,0.,0.]];load_values;[[0.,0.,0.]]] # add start and end nodes

        nodes = [Node([0.,0.,0.],Bool[0,0,0,1,1,1]); [Node([pos,0.,0.],Bool[1,1,1,1,1,1]) for pos in positions]; Node([length,0.,0.],Bool[1,1,0,1,1,1])] # singular exception can be removed by replacing with Bool[1,0,0,0,1,0]    
        loads = Asap.AbstractLoad[NodeForce(nodes[i], load_values[i]) for i in 1:lastindex(nodes)]
        elements = [Element(nodes[n-1], nodes[n], section) for n in 2:lastindex(nodes)]

    elseif is_lineload(element_loads[1])

        n_nodes = 100
        uniform_load = all_load_values[1]
        beam_load = uniform_load * length

        positions = collect(0:length/(n_nodes + 1):length) # m
        load_values = [[[0.,0.,0.]];[beam_load ./ n_nodes for i in 1:n_nodes];[[0.,0.,0.]]]

        nodes = [Node([pos,0.,0.],Bool[1,1,1,0,0,0]) for pos in positions] # singular exception can be removed by replacing with Bool[1,0,0,0,1,0]    
        nodes[1].dof = Bool[0,0,0,1,1,1]
        nodes[end].dof = Bool[1,0,0,1,1,1]
        loads = Asap.AbstractLoad[NodeForce(nodes[i], load_values[i]) for i in 1:lastindex(nodes)]
        elements = [Element(nodes[n-1],nodes[n], section) for n in 2:lastindex(nodes)]

    end

    beam_model = Asap.Model(nodes, elements, loads)

    Asap.planarize!(beam_model,:ZX)

    try
        Asap.solve!(beam_model, reprocess=true)
    catch e
        println("DEFLECTION ERROR", e)
    end

    # initial_deflection = maximum(abs.(beam_model.u[3:6:end])) # in
    # println("Initial deflection: ", initial_deflection, " in")
    param_vars = FrameVariable[]

    for element in beam_model.elements
        push!(param_vars, SectionVariable(element, section.A, 0., 1e16, :A))
        push!(param_vars, SectionVariable(element, section.Ix, 0., 1e16, :Ix))
        push!(param_vars, SectionVariable(element, section.Iy, 0., 1e16, :Iy))
        push!(param_vars, SectionVariable(element, section.J, 0., 1e16, :J))
    end

    beam_params = FrameOptParams(beam_model, param_vars);

    return beam_params

end
