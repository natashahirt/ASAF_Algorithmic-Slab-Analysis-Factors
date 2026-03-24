
"""
    plot_beam_3d(vars::Vector{T})

takes the output of the optimization (result.minimizer), turns the dimensions into sections,
and plots them.

# Arguments
- vars::Vector{T} : output variables from optimization (r.minimizer)
- beam_forces::InternalForces : the problem we were originally analyzing
- i_checks::Vector{Int} : vector of indices that the sections we're optimizing over are located at
"""
function plot_beam_3d(vars::Vector{T}, beam_forces::InternalForces, i_check::Vector{Int}) where T <: Real

    fig=Figure();
    ax=Axis3(fig[1,1],aspect=:data,title="Sections");

    beam_x = beam_forces.x

    n_vars = 4
    n_sections = Int(length(vars) / n_vars) # 4 variables to optimize over for every I_symm
    n_samples = length(beam_x)

    if i_check == []
        vars = repeat(vars,outer=2)
        n_sections = 2
        i_check = [1,n_samples]
    end

    sections = I_symm[]

    # get the areas of each section under consideration
    for i in 1:n_sections
        section = I_symm(vars[(i-1) * n_vars + 1 : (i-1) * n_vars + 4]...)
        section.x = beam_x[i_check[i]]
        push!(sections, section)
        x = [section.x for _ in section.flat_coords]
        y = [coord[1] for coord in section.flat_coords]
        z = [coord[2] for coord in section.flat_coords]
        lines!(ax,x,y,z,color=:black)
    end
    
    # draw interpolations
    for i in 1:lastindex(sections[1].flat_coords)
        interp_points = catmull_rom_interpolation([[section.x, section.flat_coords[i][1],section.flat_coords[i][2]] for section in sections], 10, closed=false)
        x = [point[1] for point in interp_points] 
        y = [point[2] for point in interp_points] 
        z = [point[3] for point in interp_points] 
        lines!(ax,x,y,z,color=:lightblue)
    end

    display(fig)

end

"""
    plot_minimizer(vars::Vector{T}, beam_forces::InternalForces, i_check::Vector{Int}) where T <: Real

plot a dashboard that shows the moment envelopes, the shear envelopes, and the geometry variables

# Arguments
- vars::Vector{T} : output variables from optimization (r.minimizer)
- beam_forces::InternalForces : the problem we were originally analyzing
- i_checks::Vector{Int} : vector of indices that the sections we're optimizing over are located at
"""

function plot_minimizer(vars::Vector{T}, beam_forces::InternalForces, i_check::Vector{Int}) where T <: Real
    
    n_vars = 4
    n_samples = length(beam_forces.x)
    n_sections = Int(length(vars)/n_vars)

    if i_check == []
        vars = repeat(vars,outer=2)
        n_sections = 2
        i_check = [1,n_samples]
    end

    fig = Figure();
    superlabel = Label(fig[0,1:2],"Summary Check")
    ax_V = Axis(fig[1:2,1],title="V");
    ax_M = Axis(fig[3:4,1],title="M");
    ax_h = Axis(fig[1,2],title="h");
    ax_w = Axis(fig[2,2],title="w");
    ax_tw = Axis(fig[3,2],title="tw");
    ax_tf = Axis(fig[4,2],title="tf");

    hlines!(ax_V, 0, color=:black)
    lines!(ax_V, beam_forces.x, beam_forces.Vy)
    scatter!(ax_V, beam_forces.x[i_check], beam_forces.Vy[i_check],color=:red)
    hlines!(ax_M, 0, color=:black)
    lines!(ax_M, beam_forces.x, beam_forces.My)
    scatter!(ax_M, beam_forces.x[i_check], beam_forces.My[i_check],color=:red)

    sections = [I_symm(vars[(i-1) * n_vars + 1 : (i-1) * n_vars + 4]...) for i in 1:n_sections]
    section_x = beam_forces.x[i_check]
    section_h = [section.h for section in sections]
    section_w = [section.w for section in sections]
    section_tw = [section.tw for section in sections]
    section_tf = [section.tf for section in sections]
    
    for i in 1:n_sections

        section = sections[i]
        h,w,tw,tf = get_geometry_vars(section)
        x = beam_forces.x[i_check[i]]

        lines!(ax_V,[x,x],[-section.Vn,section.Vn],color=:lightblue)
        lines!(ax_M,[x,x],[0,section.Mn],color=:lightblue)
        lines!(ax_h,[x,x],[0,h])
        lines!(ax_w,[x,x],[0,w])
        lines!(ax_tw,[x,x],[0,tw])
        lines!(ax_tf,[x,x],[0,tf])

    end

    res = 100

    # check if it can handle smapling random points
    sample_points = [[beam_forces.x[1]]; sort(rand(1:maximum(beam_forces.x),20)); [beam_forces.x[end]]]

    lines!(ax_h, catmull_rom_interpolation(section_x, section_h, 100, output_dimensions=2)..., color=:lightgreen)
    scatterlines!(ax_h, sample_points, [catmull_rom_point(sample_point, section_x, section_h) for sample_point in sample_points])

    lines!(ax_w, catmull_rom_interpolation(section_x, section_w, 100, output_dimensions=2)..., color=:lightgreen)
    scatterlines!(ax_w, sample_points, [catmull_rom_point(sample_point, section_x, section_w) for sample_point in sample_points])

    lines!(ax_tw, catmull_rom_interpolation(section_x, section_tw, 100, output_dimensions=2)..., color=:lightgreen)
    scatterlines!(ax_tw, sample_points, [catmull_rom_point(sample_point, section_x, section_tw) for sample_point in sample_points])

    lines!(ax_tf, catmull_rom_interpolation(section_x, section_tf, 100, output_dimensions=2)..., color=:lightgreen)
    scatterlines!(ax_tf, sample_points, [catmull_rom_point(sample_point, section_x, section_tf) for sample_point in sample_points])

    # sample single point
    sample_point = 300
    scatter!(ax_h, sample_point, catmull_rom_point(sample_point, section_x, section_h), color=:red)

    display(fig)

end
