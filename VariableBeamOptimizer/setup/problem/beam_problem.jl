"""
    beam_problem(L::Real, w_D::Real, w_L::Real; init_section::I_symm=w18x30, resolution::Int64=200, type::Symbol=:simply_supported)

set up a beam problem with typical variables, including:

# Arguments
- L::Real : length of beam [in]
- w_D::Real : dead load [kip/in]
- w_L::Real : live load [kip/in]
- init_section::I_symm=w18x30 : initial section, doesn't really matter
- resolution::Int64=200 : resolution of samples (this is for the check of the interpolation)
- type::Symbol=:simply_supported : samples from problem_types dictionary (below) to get the boundary conditions
- planar::Bool=false : should the model be planarized?
"""

function beam_problem(L::Real, w_D::Real, w_L::Real; init_section::I_symm=w18x30, resolution::Int64=200, type::Symbol=:simply_supported, planar::Bool=false)

    fixities = problem_types[type]

    w_total = 1.2w_D + 1.6w_L
    nodeStart = Node([0,0,0.], fixities[1]) # pin
    nodeEnd = Node([L,0,0.], fixities[2]) # roller
    beam = Element(nodeStart,nodeEnd,to_ASAP_section(init_section))
    load = LineLoad(beam,[0,0,-w_total]) # downward forces

    asap_model = Asap.Model([nodeStart,nodeEnd],[beam],[load])
    # planarize?
    if planar == true
        planarize!(asap_model,:XZ)
    end

    # get the moment and shear envelopes
    solve!(asap_model)

    beam_forces = InternalForces(asap_model.elements[1], asap_model, resolution=resolution);

    return beam_forces

end

const problem_types = Dict(
    :simply_supported => [:pinned,:xfree],
    :fixed => [:fixed,:fixed],
    :cantilever => [:fixed,:free]
)
