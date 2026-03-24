# Include necessary modules for slab analysis
include("utils/_utils.jl")
if get(ENV, "SLABDESIGN_SKIP_MAKIE", "0") != "1"
    include("plotting/_plotting.jl")
end
include("composite/composite.jl")
include("slab_geometries/_slab_geometries.jl")
include("slab_generation/_slab_generation.jl")
include("slab_analysis/_slab_analysis.jl")