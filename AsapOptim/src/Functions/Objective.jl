"""
    solve_truss(values::Vector{Float64}, p::TrussOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solve_truss(values::Vector{Float64}, p::TrussOptParams; linsolve_alg = UMFPACKFactorization())
    
    #populate values
    X = p.indexer.activeX ? add_values_buffer(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values_buffer(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values_buffer(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z
    A = p.indexer.activeA ? replace_values_buffer(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A


    # vₑ: 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    l = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, l)

    # Γ
    Γ = r_truss(n)

    # kₑ
    kₑ = k_truss.(p.E, A, l)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u(K, p, linsolve_alg)

    # U
    U = replace_values_buffer(zeros(p.n), p.freeids, u)

    # Store values for continuity in gradients
    return TrussResults(X,
        Y,
        Z,
        A,
        l,
        Kₑ,
        Γ,
        U)
end

"""
    solve_truss_direct(values::Vector{Float64}, p::TrussOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solve_truss_direct(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = p.indexer.activeX ? add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z
    A = p.indexer.activeA ? replace_values(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A

    # vₑ 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    l = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, l)

    # Γ
    Γ = r_truss(n)

    # kₑ
    kₑ = k_truss.(p.E, A, l)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u_direct(K, p)

    # U
    U = replace_values(zeros(p.n), p.freeids, u)

    # Store values for continuity in gradients
    return TrussResults(X,
        Y,
        Z,
        A,
        l,
        Kₑ,
        Γ,
        U)
end

function solve_truss_direct_buffer(values::Vector{Float64}, p::TrussOptParams)
    
    #populate values
    X = p.indexer.activeX ? add_values_buffer(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values_buffer(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values_buffer(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z
    A = p.indexer.activeA ? replace_values_buffer(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A

    # vₑ 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    l = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, l)

    # Γ
    Γ = r_truss(n)

    # kₑ
    kₑ = k_truss.(p.E, A, l)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u_direct(K, p)

    # U
    U = replace_values_buffer(zeros(p.n), p.freeids, u)

    # Store values for continuity in gradients
    return TrussResults(X,
        Y,
        Z,
        A,
        l,
        Kₑ,
        Γ,
        U)
end

"""
    compliance(t::TrussResults, p::TrussOptParams)

Measure of strain energy for truss structures.
"""
function compliance(t::TrussResults, p::TrussOptParams)
    dot(t.U, p.P)
end


"""
    solve_network(values::Vector{Float64}, p::NetworkOptParams)

Solve and store all relevant intermediate variables after an analysis step. This function is the basis of ALL subsequent structural analysis
"""
function solve_network(values::Vector{Float64}, p::NetworkOptParams)
    
    #populate values
    X = p.indexer.activeX ? add_values_buffer(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values_buffer(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values_buffer(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z
    q = p.indexer.activeQ ? replace_values_buffer(p.q, p.indexer.iQ, values[p.indexer.iQg] .* p.indexer.fQ) : p.q

    # fixed nodal positions
    xyz_f = [X[p.F] Y[p.F] Z[p.F]]

    # diagonal q matrix
    Q = diagm(q)

    #solve for free positions
    xyz_n = (p.Cn' * Q * p.Cn) \ (p.Pn - p.Cn' * Q * p.Cf * xyz_f)

    X2 = replace_values_buffer(X, p.N, xyz_n[:, 1])
    Y2 = replace_values_buffer(Y, p.N, xyz_n[:, 2])
    Z2 = replace_values_buffer(Z, p.N, xyz_n[:, 3])

    # Store values for continuity in gradients
    return NetworkResults(X2,
        Y2,
        Z2,
        q)
end

function target(r::NetworkResults, p::NetworkOptParams)
    norm([(p.X - r.X) (p.Y - r.Y) (p.Z - p.Z)])
end

function target(r::NetworkResults, target::Matrix{Float64})
    norm(target .- [r.X r.Y r.Z])
end

function solve_frame(values::Vector{Float64}, p::FrameOptParams; linsolve_alg = UMFPACKFactorization())

    #populate values
    X = p.indexer.activeX ? add_values_buffer(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values_buffer(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values_buffer(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z

    A = p.indexer.activeA ? replace_values_buffer(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A
    Ix = p.indexer.activeIx ? replace_values_buffer(p.Ix, p.indexer.iIx, values[p.indexer.iIxg] .* p.indexer.fIx) : p.Ix
    Iy = p.indexer.activeIy ? replace_values_buffer(p.Iy, p.indexer.iIy, values[p.indexer.iIyg] .* p.indexer.fIy) : p.Iy
    J = p.indexer.activeJ ? replace_values_buffer(p.J, p.indexer.iJ, values[p.indexer.iJg] .* p.indexer.fJ) : p.J

    # vₑ: 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    L = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, L)

    # Γ
    Γ = r_frame(n, p.Ψ)

    # kₑ
    kₑ = k_frame.(p.E, p.G, A, L, Ix, Iy, J)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u(K, p, linsolve_alg)

    # U
    U = replace_values_buffer(zeros(p.n), p.freeids, u)

    return FrameResults(
        X,
        Y,
        Z,
        A,
        Ix,
        Iy,
        J,
        L,
        Kₑ,
        Γ,
        U
    )
end

function solve_frame_direct(values::Vector{Float64}, p::FrameOptParams)

    #populate values
    X = p.indexer.activeX ? add_values(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z

    A = p.indexer.activeA ? replace_values(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A
    Ix = p.indexer.activeIx ? replace_values(p.Ix, p.indexer.iIx, values[p.indexer.iIxg] .* p.indexer.fIx) : p.Ix
    Iy = p.indexer.activeIy ? replace_values(p.Iy, p.indexer.iIy, values[p.indexer.iIyg] .* p.indexer.fIy) : p.Iy
    J = p.indexer.activeJ ? replace_values(p.J, p.indexer.iJ, values[p.indexer.iJg] .* p.indexer.fJ) : p.J

    # vₑ: 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    L = get_element_lengths(v)

    # vnₑ
    n = get_normalized_element_vectors(v, L)

    # Γ
    Γ = r_frame(n, p.Ψ)

    # kₑ
    kₑ = k_frame.(p.E, p.G, A, L, Ix, Iy, J)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u_direct(K, p)

    # U
    U = replace_values(zeros(p.n), p.freeids, u)

    return FrameResults(
        X,
        Y,
        Z,
        A,
        Ix,
        Iy,
        J,
        L,
        Kₑ,
        Γ,
        U
    )
end

function solve_frame_Pf(values::Vector{Float64}, p::FrameOptParams; linsolve_alg = UMFPACKFactorization(), dead_load::Real=0.0)

    #populate values
    X = p.indexer.activeX ? add_values_buffer(p.X, p.indexer.iX, values[p.indexer.iXg] .* p.indexer.fX) : p.X
    Y = p.indexer.activeY ? add_values_buffer(p.Y, p.indexer.iY, values[p.indexer.iYg] .* p.indexer.fY) : p.Y
    Z = p.indexer.activeZ ? add_values_buffer(p.Z, p.indexer.iZ, values[p.indexer.iZg] .* p.indexer.fZ) : p.Z

    A = p.indexer.activeA ? replace_values_buffer(p.A, p.indexer.iA, values[p.indexer.iAg] .* p.indexer.fA) : p.A
    Ix = p.indexer.activeIx ? replace_values_buffer(p.Ix, p.indexer.iIx, values[p.indexer.iIxg] .* p.indexer.fIx) : p.Ix
    Iy = p.indexer.activeIy ? replace_values_buffer(p.Iy, p.indexer.iIy, values[p.indexer.iIyg] .* p.indexer.fIy) : p.Iy
    J = p.indexer.activeJ ? replace_values_buffer(p.J, p.indexer.iJ, values[p.indexer.iJg] .* p.indexer.fJ) : p.J

    # vₑ: 
    v = get_element_vectors(X, Y, Z, p)

    # Lₑ
    L = get_element_lengths(v)
    
    # get loads
    unit_loads = A .* dead_load

    # get fixed end loads
    Pf = model_to_Pf(p, L, unit_loads)
    
    # vnₑ
    n = get_normalized_element_vectors(v, L)

    # Γ
    Γ = r_frame(n, p.Ψ)

    # kₑ
    kₑ = k_frame.(p.E, p.G, A, L, Ix, Iy, J)

    # Kₑ = ΓᵀkₑΓ
    Kₑ = get_global_ks(Γ, kₑ)

    # K
    K = assemble_global_K(Kₑ, p)

    # K⁻¹P
    u = solve_u_Pf(K, p, Pf, linsolve_alg)

    # U
    U = replace_values_buffer(zeros(p.n), p.freeids, u)

    return FrameResults(
        X,
        Y,
        Z,
        A,
        Ix,
        Iy,
        J,
        L,
        Kₑ,
        Γ,
        U
    )
end

# Function to convert unit loads to a global fixed end forces vector
function model_to_Pf(p::AbstractOptParams, lengths::Vector{Float64}, unit_loads::Vector{Float64})
    # Initialize global Pf vector with zeros
    # Assuming each element contributes 6 DOFs
    elements = p.model.elements
    total_dofs = length(p.model.nodes) * 6
    Pf = zeros(Float64, total_dofs)

    for (i, element) in enumerate(elements)
        # For each element, get the corresponding unit loads
        unit_load = unit_loads[i]

        # Calculate fixed end forces
        forces = q_local(unit_load, element)

        # Map local forces to global Pf
        nodeStart_id = element.nodeStart.nodeID
        nodeStart_globalidx = (nodeStart_id-1)*6 .+ collect(1:6)
        nodeEnd_id = element.nodeEnd.nodeID
        nodeEnd_globalidx = (nodeEnd_id-1)*6 .+ collect(1:6)

        Pf = add_values(Pf, nodeStart_globalidx, forces[1:6])
        Pf = add_values(Pf, nodeEnd_globalidx, forces[7:end])
    end

    return Pf
end

using ChainRulesCore

"""
    q_local(unit_load::Float64, element::FrameElement)

Equivalent fixed end forces for a line load given a unit load in the z-axis.
"""
function q_local(unit_load::Float64, element::Asap.FrameElement)
    LCS = element.LCS
    l = element.length

    # Load vector in LCS, assuming the load is only in the z-axis
    plocal = element.R[1:3, 1:3] * [0.0, 0.0, unit_load] .* LCS

    # Axial end forces
    ax1 = ax2 = - dot(plocal[1], LCS[1]) * l / 2

    # Perpendicular load in local y
    py = - dot(plocal[3], LCS[3])
    vy1 = vy2 = py * l / 2 # Shears in Y
    mz1 = py * l^2 / 12    # Moment 1 in Z
    mz2 = -mz1             # Moment 2 in Z

    # Perpendicular load in local z
    pz = -dot(plocal[2], LCS[2])
    vz1 = vz2 = pz * l / 2 # Shears
    my1 = -pz * l^2 / 12   # Moment 1 in Y
    my2 = -my1             # Moment 2 in Y

    return [ax1, vy1, vz1, 0.0, my1, mz1, ax2, vy2, vz2, 0.0, my2, mz2]
end

# Define the custom rule for q_local
function ChainRulesCore.rrule(::typeof(q_local), unit_load::Float64, element::Asap.FrameElement)
    y = q_local(unit_load, element)
    
    function q_local_pullback(Δy)
        # Calculate the gradient with respect to unit_load
        LCS = element.LCS
        l = element.length
        R = element.R[1:3, 1:3]
        
        # Compute the gradient of the output with respect to the unit load
        dplocal_dunit_load = R[:, 3] .* LCS
        dax1_dunit_load = -dot(dplocal_dunit_load[1], LCS[1]) * l / 2
        dpy_dunit_load = -dot(dplocal_dunit_load[2], LCS[2])
        dpz_dunit_load = -dot(dplocal_dunit_load[3], LCS[3])
        
        dvy1_dunit_load = dpy_dunit_load * l / 2
        dmz1_dunit_load = dpy_dunit_load * l^2 / 12
        dvz1_dunit_load = dpz_dunit_load * l / 2
        dmy1_dunit_load = -dpz_dunit_load * l^2 / 12
        
        # Sum the contributions from each component
        dunit_load = sum(Δy .* [dax1_dunit_load, dvy1_dunit_load, dvz1_dunit_load, 0.0, dmy1_dunit_load, dmz1_dunit_load,
                                dax1_dunit_load, dvy1_dunit_load, dvz1_dunit_load, 0.0, dmy1_dunit_load, dmz1_dunit_load])
        
        return NoTangent(), dunit_load, NoTangent()
    end
    
    return y, q_local_pullback
end