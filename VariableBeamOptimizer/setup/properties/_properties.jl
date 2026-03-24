"""
section property calculations (pure geometry)
"""

# children
include("section_geo.jl")   # geometry functions and shared properties
include("section_A.jl")     # area
include("section_Ix.jl")    # moment of inertia X axis
include("section_Iy.jl")    # moment of inertia Y axis
include("section_J.jl")     # moment of inertia polar
include("section_Sx.jl")    # elastic section modulus X axis
include("section_Sy.jl")    # elastic section modulus Y axis
include("section_Zx.jl")    # plastic seciton modulus X axis
include("section_Zy.jl")    # plastic section modulus Y axis