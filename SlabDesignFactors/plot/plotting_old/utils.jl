# clear the 0 numbers and replace with NaN
cull(x) = (x isa Number && x == 0.) ? NaN : x

# colours
ibm_colors = Dict(
    :ibm_cornflower => colorant"#648fff",
    :ibm_slateblue => colorant"#785ef0",
    :ibm_cerise => colorant"#dc267f",
    :ibm_orange => colorant"#fe6100",
    :ibm_gold => colorant"#ffb000",
    :ibm_colormap => [:ibm_cornflower, :ibm_slateblue, :ibm_cerise, :ibm_orange, :ibm_gold]
)

色 = Dict(
    :powderblue => colorant"#aeddf5",
    :skyblue => colorant"#70cbfd",
    :gold => colorant"#df7e00",
    :magenta => colorant"#dc267f",
    :orange => colorant"#e04600",
    :ceruleanblue => colorant"#00AEEF",
    :charcoalgrey => colorant"#3e3e3e",
    :irispurple => colorant"#4c2563",
    :darkpurple => colorant"#130039",
    :lilac => colorant"#A678B5", 
)

# EXTRACT THE ROTATION OF THE SYMBOL
function get_vector_1d_angle(vector_1d::Vector{<:Real}; degrees=false, clockwise=false)

    angle = atan(vector_1d[2], vector_1d[1])
 
    if clockwise == true
        if angle < 0
            angle = abs(angle)
        else
            angle = 2pi - angle
        end
    end

    if degrees == true
        angle_degrees = angle * 180/pi
        return angle_degrees
    end

    return angle

end

set_theme!(fonts = (; regular = "Arial"), titlefont = "Arial")