# When SLABDESIGN_SKIP_MAKIE=1 (see SlabDesignFactors/scripts/_scripts.jl), CairoMakie is not loaded;
# use stubs so TributaryAreas can load on headless nodes. With plotting enabled, use the Makie implementation.
if get(ENV, "SLABDESIGN_SKIP_MAKIE", "0") == "1"
    include("plotting_utils_headless.jl")
else
    include("plotting_utils_makie.jl")
end
