"""
estimation
"""

function estimate_tensile_force(Mu, d)
    return Mu/(4*d)
end

function estimate_As(Mu, d, fy)
    return estimate_tensile_force(Mu, d) / fy
end

"""
tension reinforcement design
"""

function get_β1(fc′)
    return max(0.85 - 0.05 * max(fc′ - 4, 0), 0.65) # ksi, can also be simplified to 0.85 or 0.8
end

safe_sqrt(x) = x < 0 ? 0.0 : sqrt(x)

# Using quadratic equation derived from Whitney stress block analysis
function get_Whitney_ρ(β1, R, fc′, fy)
    ρ = (1 - safe_sqrt(1 - (2/β1) * R / fc′)) / ((1/β1) * fy / fc′)
    return ρ
end

function R_design(Mn, b, d)
    R = Mn / (b * d^2)
    return R
end

function singly_reinforced_required_ρ(Mn, b, d, fy, fc′)
    β1 = get_β1(fc′)
    R = R_design(Mn, b, d)
    ρ = get_Whitney_ρ(β1, R, fc′, fy)
    return ρ
end

# steel strain at ultimate is at least 0.005; maximum allowable steel ratio to ensure ductile behavior
# this function works for ϵt = 0.005 (ACI tension limit) and ϵc = 0.004 (ACI compression limit/concrete-crushing strain)
# c/d is ratio of compression zone depth to effective depth for tension-controlled section
function singly_reinforced_max_ρ(fc′, fy; c_over_d = 0.375)
    β1 = get_β1(fc′)
    ρ = 0.85 * β1 * c_over_d * fc′/fy
    return ρ
end

function singly_reinforced_min_ρ(fc′, fy) # minimum area, ACI eqn (10 — 3)
    # if flexural strength of cracked section < uncracked section it will fail immediately upon cracking
    min_1 = (3 * sqrt(fc′ * 1e3))/(fy * 1e3) # prevent flexural failure
    min_2 = 200/(fy * 1e3) # prevent shrinkage and cracking
    return max(min_1, min_2)
end

# get moment capacity of tensile and compression steel
function get_Mu₁(Mu, Mu₂, d, d′, fy; ϕ = 0.9)
    Mu₁ = Mu - Mu₂
    As′ = Mu₁ / (ϕ * fy * (d - d′))
    return Mu₁, As′
end

# get moment capacity of tensile and compression steel
function get_As′(Mu, Mu₂, d, d′, fy; ϕ = 0.9)
    Mu₁ = Mu - Mu₂
    As′ = Mu₁ / (ϕ * fy * (d - d′))
    return As′
end

function get_Mu₂(max_ρ, b, d, fy, fc′; ϕ = 0.9) # capacity of maximum allowed tensile steel OR the capacity of tensile steel
    Mn = max_ρ * fy * b * d^2 * (1 - 0.59 * max_ρ * (fy/fc′))
    return ϕ * Mn
end

function get_longitudinal_As(b, d, d′, Mu, material::ReinforcedConcrete)

    fc′ = material.fc′
    fy = material.fy

    # find the tensile reinforcement
    ρ_single_required = singly_reinforced_required_ρ(Mu, b, d, fy, fc′)
    ρ_single_max = singly_reinforced_max_ρ(fc′, fy)
    ρ_single_min = singly_reinforced_min_ρ(fc′, fy)

    # find the compression reinforcement
    if ρ_single_max < ρ_single_required
        Mu₂ = get_Mu₂(ρ_single_max, b, d, fy, fc′)
        Mu₁, required_As′ = get_Mu₁(Mu, Mu₂, d, d′, fy)
        required_As′ = max(required_As′, 0.0)
        required_As = ρ_single_max * b * d + required_As′
    elseif ρ_single_min > ρ_single_required
        required_As′ = 0.0
        required_As = ρ_single_min * b * d
    else
        required_As′ = 0.0
        required_As = ρ_single_required * b * d
    end

    return required_As, required_As′

end

function get_longitudinal_As_differentiable(b, d, d′, Mu, material::ReinforcedConcrete)

    fc′ = material.fc′
    fy = material.fy

    # find the tensile reinforcement
    ρ_single_required = singly_reinforced_required_ρ(Mu, b, d, fy, fc′)
    ρ_single_max = singly_reinforced_max_ρ(fc′, fy)
    ρ_single_min = singly_reinforced_min_ρ(fc′, fy)

    # stable differentiable blending using sigmoids
    k = 50.0
    σ(x) = 1 / (1 + exp(-clamp(x, -50.0, 50.0)))
    s1 = σ(k * (ρ_single_required - ρ_single_max))
    s2 = σ(k * (ρ_single_min - ρ_single_required))
    w1 = s1
    w2 = (1 - s1) * s2
    w3 = (1 - s1) * (1 - s2)

    # case 1 calculations
    Mu₂_1 = get_Mu₂(ρ_single_max, b, d, fy, fc′)
    As′_1 = max(get_As′(Mu, Mu₂_1, d, d′, fy), 0.0)
    As_1 = ρ_single_max * b * d + As′_1

    # case 2 calculations  
    As′_2 = 0.0
    As_2 = ρ_single_min * b * d

    # case 3 calculations
    As′_3 = 0.0
    As_3 = ρ_single_required * b * d

    # weighted sum for smooth transition
    required_As′ = w1 * As′_1 + w2 * As′_2 + w3 * As′_3
    required_As = w1 * As_1 + w2 * As_2 + w3 * As_3

    return required_As, required_As′

end