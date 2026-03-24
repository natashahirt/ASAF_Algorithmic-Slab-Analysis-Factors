"""
    get_global_ks(rs::Vector{Matrix{Float64}}, ks::Vector{Matrix{Float64}})

Get a vector of elemental stiffness matrices in GCS given a vector of transformation matrices and a vector of elemental stiffness matrices in LCS
"""
function get_global_ks(rs::Vector{Matrix{Float64}}, ks::Vector{Matrix{Float64}})
    Symmetric.(transpose.(rs) .* ks .* rs)
end

# function get_global_ks(rs::Vector{Matrix{Float64}}, ks::Vector{Matrix{Float64}})
#     # transpose.(rs) .* ks .* rs
#     [Symmetric(r' * k * r) for (r, k) in zip(rs, ks)]
# end

"""
    assemble_global_K(elementalKs::Vector{Matrix{Float64}}, p::TrussOptParams)

Assemble the global stiffness matrix from a vector of elemental stiffness matrices
"""
function assemble_global_K(elementalKs::Vector{Symmetric{Float64,Matrix{Float64}}}, p::AbstractOptParams)

    nz = zeros(p.nnz)

    for (k, i) in zip(elementalKs, p.inzs)
        nz[i] .+= vec(k)
    end

    SparseMatrixCSC(p.n, p.n, p.cp, p.rv, nz)
end

"""
The sensitivity of K w/r/t an elemental Kₑ is the proportional stiffness added to K from Kₑ

Output is a vector: [nElements × [nDOFe × nDOFe]] of elemental stiffness matrix sensitivites
"""
function ChainRulesCore.rrule(::typeof(assemble_global_K), Eks::Vector{Symmetric{Float64,Matrix{Float64}}}, p::AbstractOptParams)
    K = assemble_global_K(Eks, p)

    function K_pullback(K̄)
        dEks = [K̄[id, id] for id in p.dofids]

        return NoTangent(), dEks, NoTangent()
    end

    return K, K_pullback
end