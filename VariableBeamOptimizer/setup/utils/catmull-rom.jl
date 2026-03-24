"""
centripetal catmull-rom interpolation

function takes a set of points and interpolates the cubic splines between them.
solves issues associated with the regular cubic spline (loops, self-intersection etc.) by
introducing a parameter α.

(after https://en.wikipedia.org/wiki/Centripetal_Catmull-Rom_spline)
"""

# helpers

"""
    get_ghost_endpoints(points::Vector{Vector{T}}) 

adds "ghost endpoints" to the list so that catmull-rom interpolates between the intended
first item and last item of the list.

# Arguments
- points::Vector{Vector{T}} : [x,y,z] points
"""
function get_ghost_endpoints(points::Vector{Vector{T}}) where T <: Real

    pushfirst!(points, 2 * points[1] - points[2])
    push!(points, 2 * points[end] - points[end-1])

    return points

end

"""
    get_quads(points::Vector{Vector{T}}; closed=false)

takes x, y, z points (or just x, y points if 2D; the script will add a z=0.) and splits them
into quadruples. The beginning and end nodes will not be included in the interpolation.

# Arguments
- positions::Vector : e.g. [[x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]]
- closed::Bool = false : if we should loop back to the beginning or not
"""
function get_quads(points::Vector{Vector{T}}; closed::Bool, i::Int=-1) where T <: Real

    quads = Vector[]

    if i > 0

        if i == 1 # ghost startpoint

            pushfirst!(points, 2 * points[1] - points[2]) 
            i += 1

        elseif i == lastindex(points)-1 # ghost endpoint

            push!(points, 2 * points[end] - points[end-1])

        end

        push!(quads, points[i-1:i+2])

    else

        if closed == false

            points = get_ghost_endpoints(points)

            for i in 1:lastindex(points)-3
                push!(quads, points[i:i+3])
            end

        else

            for i in 1:lastindex(points)
                push!(quads, [points[i],points[mod1(i+1,lastindex(points))],points[mod1(i+2,lastindex(points))],points[mod1(i+3,lastindex(points))]])
            end

        end

    end

    return quads

end

"""
    generate_t(t_last::Real, P_last::Vector{T}, P_current::Vector{T}; α::Real)

generates the centripetal knot value given the previous param and the points of interest.

# Arguments
- t_last::Real : the previous knot
- P_last::Vector{T} : the previous point
- P_current::Vector{T} : the current point
- α::Real : controls the tightness of the interpolation ∈ [0,1]
"""
function generate_t(t_last::Real, P_last::Vector{T}, P_current::Vector{T}; α::Real) where T <: Real

    euclidean_distance = sqrt(sum([(P_current[i] - P_last[i])^2 for i in 1:3]))
    t_current = (euclidean_distance)^α + t_last

    return t_current

end

# interpolator

"""
    catmull_rom_spline(quad::Vector{Vector{T}}, resolution::Int; α::Real)

gives # resolution of interpolated points between the middle two points of the quad

# Arguments
- quad::Vector{Vector{T}} : the four points (middle two are interpolated between)
- resolution::Int : how many samples we take
- α::Real : controls the tightness of the interpolation ∈ [0,1]
"""
function catmull_rom_spline(quad::Vector{Vector{T}}, resolution::Int; α::Real, output_dimensions::Int=3) where T <: Real

    if output_dimensions < 2
        interpolated_points = zeros(resolution)
    else
        interpolated_points = Vector[]
    end

    # generate t's for the quad (knot sequence t₁, t₂, t₃, t₄ from the points in the quad. t₁ always = 0.
    t₁ = 0
    t₂ = generate_t(t₁, quad[1], quad[2], α=α)
    t₃ = generate_t(t₂, quad[2], quad[3], α=α)
    t₄ = generate_t(t₃, quad[3], quad[4], α=α)

    t_range = range(t₂, t₃, length=resolution)

    for i in 1:resolution

        t = t_range[i]

        # get the intermediate function values
        A₁ = ((t₂ - t) / (t₂ - t₁)) * quad[1] + ((t - t₁) / (t₂ - t₁)) * quad[2]
        A₂ = ((t₃ - t) / (t₃ - t₂)) * quad[2] + ((t - t₂) / (t₃ - t₂)) * quad[3]
        A₃ = ((t₄ - t) / (t₄ - t₃)) * quad[3] + ((t - t₃) / (t₄ - t₃)) * quad[4]

        B₁ = ((t₃ - t) / (t₃ - t₁)) * A₁ + ((t - t₁) / (t₃ - t₁)) * A₂
        B₂ = ((t₄ - t) / (t₄ - t₂)) * A₂ + ((t - t₂) / (t₄ - t₂)) * A₃

        C = ((t₃ - t) / (t₃ - t₂)) * B₁ + ((t - t₂) / (t₃ - t₂)) * B₂

        if output_dimensions < 2
            interpolated_points[i] = C[2]
        else
            push!(interpolated_points, C)
        end

    end

    return interpolated_points

end

# full function

"""
    catmull_rom_interpolation(points::Vector{Vector{T}}; closed=false)
    (catmull_rom_interpolation(x::Vector{T}, y::Vector{T}, resolution::Int; closed=false, α=0.5, return_2d_vectors::Bool=false))

returns the catmull-rom interpolation of a series of points. Changing α changes the parameterisation
of the interpolation

(secondary method takes two values of x and y and returns the interpolated x and y vectors)

# Arguments
- positions::Vector : e.g. [[x1,y1,z1],[x2,y2,z2],...,[xn,yn,zn]]
- resolution::Int : number of points to be interpolated per spline
- closed::Bool = false : if we should loop back to the beginning or not
- α::Float :    0.5 for centripetal catmull-rom (default here)
                1 for chordal catmull-rom (curve does not follow control points tightly)
                0 for uniform catmull-rom (curve forms self-intersections and loops)
"""
function catmull_rom_interpolation(points::Vector{Vector{T}}, resolution::Int; closed=false, α=0.5, output_dimensions::Int=3, i::Int=-1) where T <: Real

    if output_dimensions == 1
        interpolated_points = zeros(resolution * (length(points)-1))
    else
        interpolated_points = Vector[]
    end

    # turn all the points into 3D
    for point in points

        @assert length(point) == 1 || length(point) == 2 || length(point) == 3 "This implementation of catmull-rom can currently only handle 1D, 2D, and 3D values."

        if length(point) == 1

            append!(point, [0,0])

        elseif length(point) == 2

            push!(point, 0)

        end

    end

    quads = get_quads(points, closed=closed, i=i)

    for i in 1:lastindex(quads)
        
        spline_points = catmull_rom_spline(quads[i], resolution, α=α, output_dimensions=output_dimensions)

        if length(quads) == 1

            interpolated_points = spline_points

        elseif output_dimensions == 1

            indices = Int.(collect(resolution * (i-1) + 1: resolution * i))
            interpolated_points[indices] = spline_points

        else

            append!(interpolated_points, spline_points)

        end

    end

    if output_dimensions == 1

        return interpolated_points

    elseif output_dimensions == 2

        x = [point[1] for point in interpolated_points]
        y = [point[2] for point in interpolated_points]

        return x, y

    else

        return interpolated_points

    end

end

function catmull_rom_interpolation(x::Vector{T}, y::Vector{T}, resolution::Int; i::Int=-1, output_dimensions::Int=3) where T <: Real

    points = [[x[i],y[i],0] for i in 1:lastindex(x)]
    output = catmull_rom_interpolation(points, resolution, closed=false, i=i, output_dimensions=output_dimensions)

    if output_dimensions == 1

        return output

    else

        x,y = output
        return x, y

    end

end

"""
    catmull_rom_point(sample::Real, x::Vector{Float64}, y::Vector{Float64})

uses find_closest to find the nearest point in the interpolation to the x you want to sample from
not particularly fast right now but workable.
enter an x variable
""" 
function catmull_rom_point(sample::Real, x::Vector{Float64}, y::Vector{Float64}; α=0.5, resolution=100)

    # verify the lengths
    if length(x) == 1

        return y[1]

    elseif length(x) == 2 # linear interpolation

        param = sample/(x[2]-x[1])
        y_sample = param * (y[2]-y[1])
        return y_sample

    end

    # if there's three points I can do a catmull-rom
    i = find_closest(sample, x)
    points = [[x[i],y[i],0] for i in 1:lastindex(x)]
    
    # generate the quad
    if i == 1 # add ghostpoint at first index

        ghostpoint = 2 * points[1] - points[2]
        quad = [ghostpoint,points[i],points[i+1],points[i+2]]

    elseif i == lastindex(x)-1 # add ghostpoint at last index

        ghostpoint = 2 * points[end] - points[end-1]
        quad = [points[i-1],points[i],points[i+1],ghostpoint]

    else

        quad = points[i-1:i+2]

    end

    # generate t's for the quad (knot sequence t₁, t₂, t₃, t₄ from the points in the quad. t₁ always = 0.
    t₁ = 0
    t₂ = generate_t(t₁, quad[1], quad[2], α=α)
    t₃ = generate_t(t₂, quad[2], quad[3], α=α)
    t₄ = generate_t(t₃, quad[3], quad[4], α=α)

    t = t₂
    t_increment = (t₃-t₂)/resolution
    
    # loop until you find the correct C
    for i in 1:resolution
        
        # get the intermediate function values
        A₁ = ((t₂ - t) / (t₂ - t₁)) * quad[1] + ((t - t₁) / (t₂ - t₁)) * quad[2]
        A₂ = ((t₃ - t) / (t₃ - t₂)) * quad[2] + ((t - t₂) / (t₃ - t₂)) * quad[3]
        A₃ = ((t₄ - t) / (t₄ - t₃)) * quad[3] + ((t - t₃) / (t₄ - t₃)) * quad[4]

        B₁ = ((t₃ - t) / (t₃ - t₁)) * A₁ + ((t - t₁) / (t₃ - t₁)) * A₂
        B₂ = ((t₄ - t) / (t₄ - t₂)) * A₂ + ((t - t₂) / (t₄ - t₂)) * A₃

        interp_point = ((t₃ - t) / (t₃ - t₂)) * B₁ + ((t - t₂) / (t₃ - t₂)) * B₂

        if interp_point[1] >= sample[1]

            return interp_point[2]

        end

        t += t_increment

    end

    sample_y = quad[3][2]
    return sample_y

end