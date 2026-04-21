function r_frame(Cxyz::AbstractArray, Ψ::Float64; tol = 1e-6)

    #local x vector cosines
    Cx, Cy, Cz = Cxyz

    #roll angle 
    cΨ = cos(Ψ)
    sΨ = sin(Ψ)

    if norm(cross(Cxyz, [0., 1., 0.])) < tol #special case for horizontal members aligned with global Y
        Λ = [0. Cy 0.;
            -Cy*cΨ 0 sΨ;
            Cy*sΨ 0 cΨ]
    else # all other
        b1 = (-Cx * Cy * cΨ - Cz * sΨ) / sqrt(Cx^2 + Cz^2)
        b2 = sqrt(Cx^2 + Cz^2) * cΨ
        b3 = (-Cy * Cz * cΨ + Cx * sΨ) / sqrt(Cx^2 + Cz^2)

        c1 = (Cx * Cy * sΨ - Cz * cΨ) / sqrt(Cx^2 + Cz^2)
        c2 = -sqrt(Cx^2 + Cz^2) * sΨ
        c3 = (Cy * Cz * sΨ + Cx * cΨ) / sqrt(Cx^2 + Cz^2)

        Λ = [Cx Cy Cz; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    [Λ zeros(3,9); zeros(3,3) Λ zeros(3,6); zeros(3,6) Λ zeros(3,3); zeros(3,9) Λ]
end

r_frame(XYZn::Matrix{Float64}, Ψ::Vector{Float64}; tol = 1e-6) = [r_frame(xyzn, psi; tol = tol) for (xyzn, psi) in zip(eachrow(XYZn), Ψ)]