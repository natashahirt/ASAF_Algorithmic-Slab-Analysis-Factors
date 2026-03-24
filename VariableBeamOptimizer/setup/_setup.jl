# abstract structures
abstract type AbstractSection end
abstract type AbstractMaterial end

# children
include("utils/_utils.jl") # math
include("materials/_materials.jl") # material properties
include("properties/_properties.jl") # functions to get section properties
include("sections/_sections.jl") # get structures for sections (e.g. I_symm)
include("problem/beam_problem.jl")
