"""
Wide-flange sections
"""

struct W_imperial 
    name::String
    A #area [in²]
    d #depth [in]
    bf #flange width [in]
    tw #web thickness [in]
    tf #flange thickness [in]
    Ix #strong moment of Inertia [in⁴]
    Zx #strong plastic modulus [in³]
    Sx #strong section modulus [in³]
    rx #strong radius of gyration [in]
    Iy #weak moment of Inertia [in⁴]
    Zy #weak plastic modulus [in³]
    Sy #weak section modulus [in³]
    ry #weak radius of gyration [in]
    J #torsional constant [in⁴]
    Cw #warping constant [in⁶]

    function W_imperial(section::W)
        
        d, bf, tw, tf, rx, ry = to_in([section.d, section.bf, section.tw, section.tf, section.rx, section.ry])
        A = to_in(section.A, power=2)
        Zx, Sx, Zy, Sy = to_in([section.Zx, section.Sx, section.Zy, section.Sy], power=3)
        Ix, Iy, J = to_in([section.Ix, section.Iy,section.J], power=4)
        Cw = to_in(section.Cw, power=6)

        depth_metric, weight_metric = [parse(Float64, x) for x in split(split(section.name, "W")[2],"X")]
        depth_imperial = Int(round(to_in(depth_metric),digits=0))
        weight_imperial = Int(round(weight_metric / 1.488,digits=0))
        if 8.5 < weight_metric / 1.488 < 9
            weight_imperial = 8.5
        end
        
        name = "W" * string(depth_imperial) * "X" * string(weight_imperial)

        return new(name, A, d, bf, tw, tf, Ix, Zx, Sx, rx, Iy, Zy, Sy, ry, J, Cw)

    end

    function W_imperial(name, A, d, bf, tw, tf, Ix, Zx, Sx, rx, Iy, Zy, Sy, ry, J, Cw)
        return new(name, A, d, bf, tw, tf, Ix, Zx, Sx, rx, Iy, Zy, Sy, ry, J, Cw)
    end
end

function to_in(x::T; power=1) where T <: Real
    return x / 25.4^power
end

function to_in(vect::Vector{T}; power=1) where T <: Real
    return [to_in(val; power=power) for val in vect]
end

function to_m(x::T; power=1) where T <: Real
    return x * 25.4^power
end

function to_m(vect::Vector{T}; power=1) where T <: Real
    return [to_m(val; power=power) for val in vect]
end

function W_imperial_to_df(self::Vector{W_imperial})
    df = DataFrame(name=String[], A=Float64[], d=Float64[], bf=Float64[], tw=Float64[], tf=Float64[], Ix=Float64[], Zx=Float64[], Sx=Float64[], rx=Float64[], Iy=Float64[], Zy=Float64[], Sy=Float64[], ry=Float64[], J=Float64[], Cw=Float64[])
    for i in 1:lastindex(self)
        @unpack name, A, d, bf, tw, tf, Ix, Zx, Sx, rx, Iy, Zy, Sy, ry, J, Cw = self[i]
        push!(df, [name, A, d, bf, tw, tf, Ix, Zx, Sx, rx, Iy, Zy, Sy, ry, J, Cw])
    end
    return df
end

function get_geometry_vars(self::W_imperial)

    return [self.d, self.bf, self.tw, self.tf]
end

W_metric_names = AsapToolkit.names[AsapToolkit.Wrange]
W_imperial_names = String[]

function W_imperial(name::String)

    if name in W_metric_names
        return W_imperial(W(name))
    elseif name in W_imperial_names
        mm_name = W_dict_imperial[name]
        return W_imperial(W(mm_name))
    end

end

const SECTIONS = [W_imperial(name) for name in AsapToolkit.names[AsapToolkit.Wrange]]
W_imperial_names = [Section.name for Section in SECTIONS]

W_dict_metric = Dict()
W_dict_imperial = Dict()

for i in 1:lastindex(W_imperial_names)
    W_dict_metric[W_metric_names[i]] = W_imperial_names[i]
    W_dict_imperial[W_imperial_names[i]] = W_metric_names[i]
end

p = sortperm([Section.A for Section in SECTIONS])

function toASAPframe_W(section::Union{W, W_imperial}, E::Real, G::Real; convert::Bool=false)
    if convert
        A = to_m(section.A, power=2)
        Ix = to_m(section.Ix, power=4)
        Iy = to_m(section.Iy, power=4)
        J = to_m(section.J, power=4)
        return Section(A, E, G, Ix, Iy, J)
    else
        return Section(section.A, E, G, section.Ix, section.Iy, section.J)
    end
end

const allW_imperial() = SECTIONS[p]