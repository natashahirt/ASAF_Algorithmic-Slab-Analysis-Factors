"""
    get_beam_releases(nodes::Vector{Node}, connectivities::Vector{Vector{Int64}}; element_ids::Vector{Int64}=Int64[])

Determines the release conditions for beams based on their connectivity and optional element IDs.

# Arguments
- `nodes`: A vector of `Node` objects.
- `connectivities`: A vector of vectors, each containing two integers representing node indices.
- `element_ids`: An optional vector of element IDs.

# Returns
- A tuple containing a vector of release symbols and a vector of new connectivities.
"""
function get_beam_releases(nodes::Vector{Node}, connectivities::Vector{Vector{Int64}}; element_ids::Vector{Int64}=Int64[], i_holes::Vector{Int64}=Int64[])
    releases = Symbol[]
    new_connectivities = Vector[]

    if !isempty(element_ids) && length(connectivities) == length(element_ids)
        connectivities_grouped = Vector[]
        
        # Group connectivities by element ID
        for i in 1:maximum(element_ids)
            for j in 1:lastindex(connectivities)
                if element_ids[j] == i
                    if length(connectivities_grouped) < i
                        push!(connectivities_grouped, [connectivities[j]])
                    else
                        push!(connectivities_grouped[i], connectivities[j])
                    end
                end
            end
        end

        for (i, group) in enumerate(connectivities_grouped)
            # Check if any nodes in the group are in i_holes
            any_holes = false
            # Count how many nodes in the group are in holes
            hole_count = 0
            for connectivity in group
                if connectivity[1] in i_holes
                    hole_count += 1
                end
                if connectivity[2] in i_holes
                    hole_count += 1
                end
            end

            any_holes = hole_count > 1

            # If any nodes are in holes, make entire group fixed-fixed
            """if any_holes
                for connectivity in group
                    push!(releases, :fixedfixed)
                    push!(new_connectivities, connectivity)
                end
                continue
            end"""

            if length(group) == 1
                release = :freefree
                connectivity = group[1]
                push!(releases, release)
                push!(new_connectivities, connectivity)
                continue
            end

            group_nodes = Int64[]
            for connectivity in group
                append!(group_nodes, connectivity)
            end

            node_counts = Dict()
            for node in group_nodes
                node_counts[node] = get(node_counts, node, 0) + 1
            end
            
            for (j, connectivity) in enumerate(group)
                nodeStart = connectivity[1]
                nodeEnd = connectivity[2]

                if node_counts[nodeStart] == 2 && node_counts[nodeEnd] == 2 
                    release = :fixedfixed
                elseif node_counts[nodeStart] == 2 && node_counts[nodeEnd] == 1
                    release = :fixedfree
                elseif node_counts[nodeStart] == 1 && node_counts[nodeEnd] == 2
                    release = :freefixed
                else
                    release = :freefree
                end

                push!(releases, release)
                push!(new_connectivities, connectivity)
            end
        end
    else
        for connectivity in connectivities
            nodeStart = nodes[connectivity[1]]
            nodeEnd = nodes[connectivity[2]]
    
            if (nodeStart.id in [:column, :wall]) && (nodeEnd.id in [:column, :wall]) 
                release = :freefree
            elseif nodeStart.id in [:column, :wall]
                release = :freefixed
            elseif nodeEnd.id in [:column, :wall]
                release = :fixedfree
            else
                release = :fixedfixed
            end
    
            push!(releases, release)
            push!(new_connectivities, connectivity)
        end
    end

    return releases, new_connectivities
end

"""
    connect_nodes(node_list::Vector{Node}, connection::Vector{Int64}, release::Symbol; section=toASAPframe(W(W_dict_imperial["W8X35"]), 210000, 81000, unit = :m))

Connects two nodes and creates an element with a specified release condition and section.

# Arguments
- `node_list`: A vector of `Node` objects.
- `connection`: A vector containing two integers representing node indices.
- `release`: A symbol representing the release condition.
- `section`: An optional section parameter.

# Returns
- An `Element` object representing the connection.
"""
function connect_nodes(node_list::Vector{Node}, connection::Vector{Int64}; release::Symbol=:fixedfixed, section::Asap.Section)
    nodeStart = node_list[connection[1]]
    nodeEnd = node_list[connection[2]]

    posStart = nodeStart.position
    posEnd = nodeEnd.position

    element = Element(nodeStart, nodeEnd, section, release=release)
    element.length = sqrt((posEnd[1]-posStart[1])^2 + (posEnd[2]-posStart[2])^2 + (posEnd[3]-posStart[3])^2)
    element.nodeStart.nodeID = connection[1]
    element.nodeEnd.nodeID = connection[2]

    return element
end

"""
    connect_nodes(node_list::Vector{Node}, connectivity_list, releases::Vector{Symbol}; sections::Vector=[])

Connects multiple nodes based on a list of connectivities and release conditions.

# Arguments
- `node_list`: A vector of `Node` objects.
- `connectivity_list`: A list of connectivities, each containing two node indices.
- `releases`: A vector of release symbols.
- `sections`: An optional vector of sections.

# Returns
- A vector of `Element` objects.
"""
function connect_nodes(node_list::Vector{Node}, connectivity_list::Union{Vector{Vector}, Vector{Vector{Int64}}}; releases::Vector{Symbol}=Symbol[], sections::Union{Vector{Any}, Vector{String}, Vector{Section}, Vector{Vector}}=[])
    E = steel_mm.E
    G = steel_mm.G
    unitfactors = [1e-6, 1, 1, 1e-12, 1e-12, 1e-12]

    if isempty(releases)
        releases = [:fixedfixed for i in 1:lastindex(connectivity_list)]
    end

    if isempty(sections)
        default_section = W(W_dict_imperial["W8X35"])
        sections = [toASAPframe(default_section, E, G, unit = :m) for i in 1:lastindex(connectivity_list)]
    elseif typeof(sections) == Vector{String}
        sections = [toASAPframe(W(W_dict_imperial[section]), E, G, unit = :m) for section in sections]
    elseif typeof(sections) == Vector{Vector}
        A, Ix, Iy, J = Float64[], Float64[], Float64[], Float64[]
        for i in 1:lastindex(sections)
            I_imperial = I_symm(sections[i]...)
            push!(A, I_imperial.A * 25.4^2) # in² -> mm²
            push!(Ix, I_imperial.Ix * 25.4^4) # in⁴ -> mm⁴
            push!(Iy, I_imperial.Iy * 25.4^4) # in⁴ -> mm⁴
            push!(J, I_imperial.J * 25.4^4) # in⁴ -> mm⁴
        end
        sections = [Section([A[i], E, G, Ix[i], Iy[i], J[i]] .* unitfactors...) for i in 1:lastindex(sections)]
    elseif typeof(sections) == Vector{Section}
        sections = sections
    end

    elements = Element[]
    for (i, connection) in enumerate(connectivity_list)
        release = releases[i]
        section = sections[i]
        if isa(connection, Set)
            connection = collect(connection)
        end
        element = connect_nodes(node_list, connection, release=release, section=section)
        push!(elements, element)
    end

    return elements
end

"""
    check_overlap(nodes::Vector{Node}, edge::Vector{Int64})

Checks for overlapping nodes along a given edge.

# Arguments
- `nodes`: A vector of `Node` objects.
- `edge`: A vector containing two integers representing node indices.

# Returns
- A vector of indices of overlapping nodes.
"""
function check_overlap(nodes::Vector{Node}, edge::Vector{Int64})
    m, c = line_equation(nodes[edge[1]], nodes[edge[2]])
    overlapping_nodes = Int64[]

    for i in 1:lastindex(nodes)
        if i in edge
            continue
        end

        if m == 0 && nodes[i].position[2] == nodes[edge[1]].position[2] # m is zero, so check for same y
            push!(overlapping_nodes, i)
        elseif isinf(m) && nodes[i].position[1] == nodes[edge[1]].position[1] # m is infinite, so check for same x
            push!(overlapping_nodes, i)
        elseif nodes[i].position[2] == m * nodes[i].position[1] + c 
            push!(overlapping_nodes, i)
        end
    end

    return overlapping_nodes
end

"""
    get_stiffness_ratio(element_1::Element, element_2::Element)

Calculates the stiffness ratio between two elements.

# Arguments
- `element_1`: The first `Element` object.
- `element_2`: The second `Element` object.

# Returns
- The stiffness ratio as a float.
"""
function get_stiffness_ratio(element_1::Element, element_2::Element)
    stiffness_1 = element_1.section.E * element_1.section.Ix / get_length(element_1)^3
    stiffness_2 = element_2.section.E * element_2.section.Ix / get_length(element_2)^3

    ratio = stiffness_1 / (stiffness_1 + stiffness_2)
    return ratio
end