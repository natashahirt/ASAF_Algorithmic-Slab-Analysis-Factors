"""
    connect_columns(nodes::Vector{Node}; i_columns=Int64[], frame_height=1.0)

Connects columns to the given nodes. If no column indices are provided, it finds the convex hull
and extracts the corner nodes as column supports. It then builds the columns and labels the elements.

# Arguments
- `nodes::Vector{Node}`: A vector of nodes to connect columns to.
- `i_columns=Int64[]`: Optional indices of nodes to be used as columns.
- `frame_height=1.0`: The height of the frame for column placement.

# Returns
- `nodes`: Updated vector of nodes with new ground nodes added.
- `column_elements`: Vector of column elements created.
"""
function connect_columns(nodes::Vector{Node}; i_columns=Int64[], frame_height=1.0, section::Asap.Section)

    # If no column indices are provided, find the convex hull and extract corner nodes
    if isempty(i_columns)
        i_hull, i_columns = quickhull_algorithm(nodes)
    end

    # Remove duplicate column indices
    i_columns = collect(Set(i_columns))

    # Initialize an empty vector to store column elements
    column_elements = Element[]

    # Iterate over each column index
    for i in i_columns
        # Create a new ground node vertically below the current column node
        ground_node = Node([nodes[i].position[1], nodes[i].position[2], nodes[i].position[3] - frame_height], :fixed)
        
        # Add the ground node to the nodes vector
        push!(nodes, ground_node)
        
        # Connect the current node to the new ground node
        column_element = connect_nodes(nodes, [i, lastindex(nodes)], release=:freefixed, section=section)

        # Set the node IDs for the column element
        column_element.nodeStart.id = :column
        column_element.nodeEnd.id = :fixed

        # Add the column element to the column_elements vector
        push!(column_elements, column_element)
    end

    return nodes, column_elements
end

"""
    set_material(element::Element, material::Material, units=:m)

Sets the material properties for an element using a material object.

# Arguments
- `element::Element`: The element to set material properties for.
- `material::Material`: The material object containing properties.
- `units=:m`: The units for the material properties (default is meters).

# Returns
- `element`: The element with updated material properties.
"""
function set_material(element::Element, material::Material, units=:m)
    E, G, ρ = material.E, material.G, material.ρ
    return set_material(element, E, G, ρ)
end

"""
    set_material(element::Element, E::Float64, G::Float64, ρ::Float64)

Sets the material properties for an element using individual property values.

# Arguments
- `element::Element`: The element to set material properties for.
- `E::Float64`: Young's modulus.
- `G::Float64`: Shear modulus.
- `ρ::Float64`: Density.

# Returns
- `element`: The element with updated material properties.
"""
function set_material(element::Element, E::Float64, G::Float64, ρ::Float64)
    element.section.E, element.section.G, element.section.ρ = E, G, ρ
    return element
end
