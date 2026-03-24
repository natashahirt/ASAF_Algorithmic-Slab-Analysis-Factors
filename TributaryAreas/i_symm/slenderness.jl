# AISC 360 Table B4.1b — Slenderness Limits
# These re-export / alias the functions already on I_symm for callers
# that prefer the (section, material) calling convention.

# get_slenderness(s::I_symm) and get_compression_factors(s::I_symm) are
# defined in the main I_symm.jl and dispatch on I_symm directly.
# No additional code needed here — this file is kept so the include
# structure in _i_symm.jl remains valid.
