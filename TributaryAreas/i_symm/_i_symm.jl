# ==============================================================================
# AISC 360 Steel Design Checks — Extensions for I_symm
# ==============================================================================
# These files add standalone design check functions that accept I_symm
# (defined in VariableBeamOptimizer) and Metal. The core struct already
# computes capacities on construction; these are for on-demand checks
# with different parameters (e.g. different Lb, Cb, axis, effective length).

include("slenderness.jl")  # Table B4.1b / E7 — slender element reduction
include("flexure.jl")      # Chapter F  — standalone Mn with custom Lb/Cb
include("shear.jl")        # Chapter G  — Cv1, strong/weak axis Vn
include("compression.jl")  # Chapter E  — Fe, Fcr, Pn
include("torsion.jl")      # Design Guide 9 — torsion for W-shapes
