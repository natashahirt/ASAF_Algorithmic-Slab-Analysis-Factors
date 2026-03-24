include("Utilities.jl")
export replace_values
export add_values

include("Geometry.jl")

include("Rtruss.jl")
include("Rframe.jl")

include("Ktruss.jl")
include("Kframe.jl")

include("K.jl")

include("Solve.jl")

include("Objective.jl")
export solve_truss
export solve_truss_direct
export compliance

export solve_network
export target

export solve_frame
export solve_frame_direct

include("PostProcessing.jl")
export axial_force
export axial_stress
export updatemodel
export updatenetwork