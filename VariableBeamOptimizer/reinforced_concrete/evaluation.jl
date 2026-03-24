# down the line for verification

function check_compression_steel_yielded(b, d, d′, As, As′, fy, fc′)
    # Check if compression steel has yielded using ACI 318-19
    # Left side of equation: (As - As′)/(b*d)
    # Right side: 0.85 * β₁ * (fc′/fy) * (87000/(87000-fy)) * (d′/d)
    # We make assumption that the compression steel yields at failure
    β1 = get_β1(fc′)
    left_side = (As - As′)/(b*d)
    right_side = 0.85 * β1 * (fc′/fy) * (87000/(87000-fy)) * (d′/d)
    
    return left_side - right_side 
end


function check_tension_controlled_section(b, d, d′, As, As′, fy, fc′)
    β1 = get_β1(fc′)
    # Check if section is tension-controlled per ACI 318-19
    c = get_c(b, d, d′, As, As′, fy, fc′, β1)
    if c <= 0 return -0.375 end
    c_over_d = c / d
    println("c/d: ", c_over_d)
    return c_over_d - 0.375
end