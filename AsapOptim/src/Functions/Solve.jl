"""
    solve_u(K::SparseMatrixCSC{Float64, Int64}, p::TrussOptParams)

Displacement of free DOFs
"""
function solve_u(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams, alg::LinearSolve.SciMLLinearSolveAlgorithm)
    linearproblem = LinearProblem(K[p.freeids, p.freeids], p.P[p.freeids])
    sol = LinearSolve.solve(linearproblem, alg)

    return sol.u
end

function ChainRulesCore.rrule(::typeof(solve_u), K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams, alg::LinearSolve.SciMLLinearSolveAlgorithm)
    linearproblem = LinearProblem(K[p.freeids, p.freeids], p.P[p.freeids])
    linsolve = LinearSolve.init(linearproblem, alg)
    sol1 = LinearSolve.solve!(linsolve)
    u = sol1.u

    function solve_u_pullback(ū)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        # dKΔ = cg(K[p.freeids, p.freeids], ū)
        linsolve.b = ū
        sol2 = LinearSolve.solve!(linsolve)

        #assign to proper indices
        @views dudK[p.freeids, p.freeids] = kron(u', sol2.u)

        return NoTangent(), -dudK, NoTangent(), NoTangent()
    end

    return u, solve_u_pullback
end

"""
    solve_u_direct(K::SparseMatrixCSC{Float64, Int64}, p::TrussOptParams)

Displacement of free DOFs using direct solve
"""
function solve_u_direct(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams)
    id = p.freeids
    K[id, id] \ p.P[id]
end

function ChainRulesCore.rrule(::typeof(solve_u_direct), K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams)
    u = solve_u_direct(K, p)

    function solve_u_direct_pullback(ū)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        dKΔ = K[p.freeids, p.freeids] \ ū

        #assign to proper indices
        dudK[p.freeids, p.freeids] .= kron(u', dKΔ)

        return NoTangent(), -dudK, NoTangent()
    end

    return u, solve_u_direct_pullback
end

"""
    solve_u_Pf(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams, Pf::Vector{Float64}, alg::LinearSolve.SciMLLinearSolveAlgorithm)

Displacement of free DOFs including fixed end forces. Solves the system K[freeids, freeids] * u = P[freeids] + Pf[freeids]
"""
function solve_u_Pf(K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams, Pf::Vector{Float64}, alg::LinearSolve.SciMLLinearSolveAlgorithm)
    linearproblem = LinearProblem(K[p.freeids, p.freeids], p.P[p.freeids] + Pf[p.freeids])
    sol = LinearSolve.solve(linearproblem, alg)

    return sol.u
end

function ChainRulesCore.rrule(::typeof(solve_u_Pf), K::SparseMatrixCSC{Float64, Int64}, p::AbstractOptParams, Pf::Vector{Float64}, alg::LinearSolve.SciMLLinearSolveAlgorithm)
    linearproblem = LinearProblem(K[p.freeids, p.freeids], p.P[p.freeids] + Pf[p.freeids])
    linsolve = LinearSolve.init(linearproblem, alg)
    sol1 = LinearSolve.solve!(linsolve)
    u = sol1.u

    function solve_u_Pf_pullback(ū)
        # Initialize
        dudK = zeros(p.n, p.n)
        dudPf = zeros(p.n)

        # Solve adjoint system for sensitivities
        linsolve.b = ū
        sol2 = LinearSolve.solve!(linsolve)
        
        # Assign sensitivities to proper indices
        @views dudK[p.freeids, p.freeids] = kron(u', sol2.u)
        @views dudPf[p.freeids] = sol2.u

        return NoTangent(), -dudK, NoTangent(), dudPf, NoTangent()
    end

    return u, solve_u_Pf_pullback
end