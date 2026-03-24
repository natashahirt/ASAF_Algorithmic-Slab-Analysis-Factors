"""
    generate_objective(section::I_symm, constant::Vector{Bool})

2 things:
- initial variables: 
    - input: h w tw tf + bool fixed or not
    - output: h w tw tf w tw tf w tw tf w tw tf only repeat the variables that are changing
    -
- objective variables: 
    - input: h w tw tf w tw tf w tw tf w tw tf
    - output: h w tw tf h w tw tf h w tw tf h w tw tf_indices
I need an index map, and I can generate that outside of the broadcast function itself (keep the broadcast differentiable)
"""
function generate_broadcast(fix_vars::Vector{Bool}, n_sections::Int64)
    # example values are for (fix_vars = [true, true, false, false], n_sections = 5)

    broadcast_indices = collect(1:n_sections * length(fix_vars)) # 1, 2, 3, 4, 5, ... , 20
    
    fixed_i = findall(x->x==true,fix_vars)
    free_i = findall(x->x==false,fix_vars)
    
    for i in fixed_i

        fix_indices = collect(0:n_sections-1) .* length(fix_vars) .+ i # get every fourth variable, for instance
        broadcast_indices[fix_indices] .= i # apply the broadcast index 

    end

    free_indices = Int[]

    for i in free_i

        append!(free_indices, collect(0:n_sections-1) .* length(fix_vars) .+ i)

    end

    free_indices = sort(free_indices)
    count = 0

    for i in free_indices
    
        count += 1

        while count in fixed_i
            
            count += 1

        end
           
        broadcast_indices[i] = count

    end
    
    function broadcast_in(vars::Vector{T}) where T <: Real # takes a 4 element vector and turns it into the design vector

        fixed_values = vars # h, w 
        repeat_values = repeat([vars[i] for i in 1:lastindex(fix_vars) if (fix_vars[i] == false)], outer=n_sections-1)

        return ([fixed_values; repeat_values])

    end

    function broadcast_out(vars::Vector{T}) where T <: Real # takes a design vector and turns it into the n_sections * 4 element vector

        broadcasted_vars = vars[broadcast_indices]
        
        return broadcasted_vars

    end
    
    return Function[broadcast_in, broadcast_out]

end


"""
    objective_I_symm(vars::Vector{Real})

objective can be :V
"""

# new
function generate_objective(params::SlabSizingParams, beam_params::FrameOptParams, beam_forces::InternalForces, broadcast_out::Function, i_check::Vector{Int}, n_sections::Int64; interpolation::Function=akima_interpolation)
    
    function I_symm_prismatic(vars::Vector{T}) where T <: Real 

        @assert length(vars) == 4 "Your prismatic design vector isn't prismatic."

        A = A_I_symm(vars...) # inches²
        
        beam_length = beam_forces.x[end] # already in inches

        # volume
        V = A * beam_length # inches³

        print(".")

        # deflection

        if params.deflection_limit == true

            δ_max = beam_length/params.serviceability_lim # generally beam_length/360, beam_length/200 if roof beam
            δ_local = get_element_deflection(beam_params, vars, material=steel_ksi)
            δ = maximum(abs.(δ_local)) # inches

            if δ - δ_max <= 0
                δ_penalty = (δ - δ_max)^2
            else
                δ_penalty = (δ - δ_max)^4
            end
            
            # print statement
            print(".")

        else
        
            δ_penalty = 0
            
            # print statement
            print(".")

        end

        # objective
        objective = V + δ_penalty

        return objective

    end

    return I_symm_prismatic

end

"""
    constraints_model(Mu::Real, Vu::Real)

generate a list of constraints for the model. Inequality are of the form f(x) <= 0
the output of f(x) can be a vector, in which case f(x) .<= 0 (elementwise). This allows us to
condense multiple geometric constraints into a single constraint function.

# Arguments
- Mu::Real : required moment
- Vu::Real : required shear
- n_sections::Int=5 : how many sections we want to optimize over by default
- i_check::Vector{Int} : the indices of the sampled sections that represent the sections we want to optimize
"""

function inequality_constraints_I_symm(beam_forces::InternalForces, n_sections::Int; i_check::Vector{Int}=[], n_validation::Int=0, interpolation_function::Function=akima_interpolation, broadcast_out::Function)

    Mu, Vu, x = beam_forces.My, beam_forces.Vy, beam_forces.x
    Vu = abs.(Vu) # absolute value
    Mu = abs.(Mu) # absolute value of moment
    
    n_vars = 4 # how many variables are in the I_symm beam
    
    # using LRFD, Rᵤ (required strength) <= ϕ * Rₙ (design strength)
    ϕ_b = 0.9 # resistance factor for moment, Mu ≤ ϕ_b * Mn
    ϕ_v = 0.9 # resistance factor for shear, Vu ≤ ϕ_v * Vn

    """ if there are no specific sections to check, size for the moment envelope"""

    if length(i_check) == 0

        functions = Function[]

        i_check = [1] 

        Mu_section = maximum(Mu)
        Vu_section = maximum(Vu) 

        function single_section(vars::Vector)

            h, w, tw, tf = vars

            constraint_vector = get_constraint_vector(h, w, tw, tf, Mu_section, Vu_section, ϕ_b=ϕ_b, ϕ_v=ϕ_v)

            return constraint_vector

        end

        push!(functions, single_section)

        return functions

    end

    """ constraints for multiple sections
    # includes both the sections we're optimizing as well
    # as series of checks (validation sections)

    functions = Function[]

    # optimization sections
    
    for i in 1:n_sections

        function valid_section(vars::Vector{T}) where T <: Real

            vars = broadcast_out(vars)

            h, w, tw, tf = vars[(i-1) * n_vars + 1 : (i-1) * n_vars + 4]

            Mu_section = Mu[i_check[i]]
            Vu_section = Vu[i_check[i]]
            
            constraint_vector = get_constraint_vector(h, w, tw, tf, Mu_section, Vu_section, ϕ_b=ϕ_b, ϕ_v=ϕ_v)

            println("constraint_vector: $constraint_vector")q

            return constraint_vector

        end

        push!(functions, valid_section)

    end

    # validation sections

    if n_validation > 0

        if n_validation == 1

            n_validation = 2

        end

        x_restriction = collect(range(minimum(x[i_check]),maximum(x[i_check]),n_validation))

        h_indices = Vector{Int}[(collect(range(1,length(i_check)*4-3, length(i_check))))][1]
        w_indices = h_indices .+ 1
        tw_indices = h_indices .+ 2
        tf_indices = h_indices .+ 3

        for x_sample in x_restriction

            function valid_section(vars::Vector{T}) where T <: Real

                vars = broadcast_out(vars)
                
                section_x = x[i_check]

                section_h = vars[h_indices]
                section_w = vars[w_indices]
                section_tw = vars[tw_indices]
                section_tf = vars[tf_indices]

                h = interpolation_function(x_sample, section_x, section_h)
                w = interpolation_function(x_sample, section_x, section_w)
                tw = interpolation_function(x_sample, section_x, section_tw)
                tf = interpolation_function(x_sample, section_x, section_tf)
                
                Mu_section = interpolation_function(x_sample, x, Mu)
                Vu_section = interpolation_function(x_sample, x, Vu)

                constraint_vector = get_constraint_vector(h, w, tw, tf, Mu_section, Vu_section, ϕ_b=ϕ_b, ϕ_v=ϕ_v)
                return constraint_vector

            end

            push!(functions, valid_section)

        end

    end

    return functions"""

end
