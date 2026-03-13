# Get data from old figure's first plot
function copy_axis(old_ax::Axis, new_ax::Axis; alpha=1.0, scatter=true, lines=true, beams=true)
    if !isempty(old_ax.scene.plots)
        for plot in old_ax.scene.plots
            if isa(plot, Scatter) && scatter
                x_y_data = plot.args
                x, y = x_y_data[1][], x_y_data[2][]
                attributes = plot.attributes  # Get all attributes
                color = attributes[:color]
                strokecolor = attributes[:strokecolor]
                strokewidth = attributes[:strokewidth]
                marker = attributes[:marker]
                markersize = attributes[:markersize]
                if !isa(color[], Symbol)
                    colormap = attributes[:colormap]
                    colorrange = attributes[:colorrange]
                    scatter!(new_ax, x, y, color=color, colormap=colormap, colorrange=colorrange, marker=marker, markersize=markersize, strokecolor=strokecolor, strokewidth=strokewidth, alpha=alpha)
                else
                    scatter!(new_ax, x, y, color=color, marker=marker, markersize=markersize, strokecolor=strokecolor, strokewidth=strokewidth, alpha=alpha)
                end
            elseif isa(plot, Lines) && lines
                x_y_data = plot.args
                x, y = x_y_data[1][], x_y_data[2][]
                attributes = plot.attributes  # Get all attributes
                color = attributes[:color]
                if !isa(color[], Symbol)
                    colormap = attributes[:colormap]
                    colorrange = attributes[:colorrange]
                    lines!(new_ax, x, y, color=color, colormap=colormap, colorrange=colorrange, alpha=alpha)
                elseif beams == true
                    lines!(new_ax, x, y, color=color, alpha=alpha)
                end
            end
        end
    end

    return new_ax
end

function get_beams(old_ax::Axis)
    beams = []
    for plot in old_ax.scene.plots
        if isa(plot, Lines) && isa(plot.attributes[:color][], Symbol)
            push!(beams, plot)
        end
    end
    return beams
end

function parse_sections(sections_str::String)
    parsed_sections = Meta.parse(sections_str)
    parsed_array = eval(parsed_sections)
    return parsed_array
end

function parse_sections(sections_vector::Vector{String})
    parsed_sections = Vector{Any}(undef, length(sections_vector))
    for (i, section) in enumerate(sections_vector)
        parsed_sections[i] = parse_sections(section)
    end
    return parsed_sections
end

function parse_ids(ids_vector::Vector{SubString{String}})
    parsed_ids = Vector{Any}(undef, length(ids_vector))
    for (i, id) in enumerate(ids_vector)
        parsed_ids[i] = replace(replace(replace(id, r"[\[\]]" => ""), "\"" => ""), "Any" => "")
        parsed_ids[i] = strip(parsed_ids[i])
    end
    return string.(parsed_ids)
end

function plot_slab(self::SlabAnalysisParams, sections::Union{Vector{String}, Vector{Float64}, Vector{Any}, Vector{Vector{Float64}}}; section_names::Vector{String}=String[], text=true, mini=false, background=true, collinear=false)

    model = self.model

        # Make a new figure and copy the axis
    if mini
        new_fig = Figure(size=(380,250))
    else
        new_fig = Figure(size=(600,400))
    end
    
    new_ax = Axis(new_fig[1,1], aspect=DataAspect())

    if background
        old_ax = self.plot_context.ax
        new_ax = copy_axis(old_ax, new_ax, alpha=0.5, scatter=false, beams=false)
    end

    hidespines!(new_ax)
    hidedecorations!(new_ax)
    
    # Copy relevant axis properties
    elements = self.model.elements[:beam]
    text_plots = []

    if typeof(sections) <: Vector{Any}
        sections = string.(sections)
        if !contains(sections[1], "W")
            sections = parse_sections(sections)
        end
    end

    if typeof(sections) <: Vector{String}
        areas = [W_imperial(section).A for section in sections]
    elseif typeof(sections) <: Vector{Vector{Float64}}
        areas = [I_symm(section...).A for section in sections]
        sections = string.(round.(areas, digits=2))
    else
        areas = sections
        sections = string.(sections)
    end

    area_range = (0, sqrt(maximum(areas)))

    # Check if area_range is valid
    if area_range[1] == area_range[2]
        # Adjust the range slightly or set a default color
        area_range = (area_range[1], area_range[1] + 1e-5)
    end

    # Plot each beam element
    for (i, element) in enumerate(elements)
        if element.nodeStart.id == :wall && element.nodeEnd.id == :wall
            lines!(new_ax, [element.nodeStart.position[1], element.nodeEnd.position[1]], [element.nodeStart.position[2], element.nodeEnd.position[2]], color=:grey, linewidth=5)
            continue
        end
        
        if areas[i] > 0

            x = [element.nodeStart.position[1], element.nodeEnd.position[1]]
            y = [element.nodeStart.position[2], element.nodeEnd.position[2]]
            linewidth = mini ? sqrt(areas[i]) / 2 : sqrt(areas[i])
            color = sqrt(areas[i])

            # Calculate clipped line coordinates by moving inward from endpoints
            # Get release type from element and determine DOFs
            release_type = typeof(element).parameters[1]

            start_dof = if release_type <: Union{Asap.FixedFixed, Asap.FixedFree}  && sum(element.nodeStart.dof[4:6]) == 0
                [0,0,0,0,0,0] # Fixed start
            else
                [0,0,0,1,1,1] # Free start
            end
            end_dof = if release_type <: Union{Asap.FixedFixed, Asap.FreeFixed}  && sum(element.nodeEnd.dof[4:6]) == 0
                [0,0,0,0,0,0] # Fixed end
            else
                [0,0,0,1,1,1] # Free end
            end

            # Default to no clipping for moment connections
            if mini
                clip_amount = .5
            else
                clip_amount = .5
            end

            clip_start = sum(start_dof[4:6]) == 0 ? 0 : clip_amount  # Clip if rotation DOF is fixed
            clip_end = sum(end_dof[4:6]) == 0 ? 0 : clip_amount      # Clip if rotation DOF is fixed
            
            dx = x[2] - x[1]
            dy = y[2] - y[1]
            element_length = sqrt(dx^2 + dy^2)
            
            # Calculate parametric distances based on connection types
            t1 = clip_start/element_length
            t2 = (element_length-clip_end)/element_length
            
            x_clipped = [x[1] + t1*dx, x[1] + t2*dx]
            y_clipped = [y[1] + t1*dy, y[1] + t2*dy]
            lines!(new_ax, x_clipped, y_clipped, linewidth=linewidth, color=color, colorrange=area_range, colormap=:BuPu)

            # Calculate the length of the beam
            element_length = sqrt(dx^2 + dy^2)

            # Text
            # Calculate midpoint coordinates
            mid_x = (x[1] + x[2]) / 2
            mid_y = (y[1] + y[2]) / 2
            
            # Calculate rotation angle in radians
            rotation = let θ = atan((y[2] - y[1]), (x[2] - x[1]))
                if abs(θ - π/4) < 0.1
                    π/4  # Positive 45 degrees
                elseif abs(θ + π/4) < 0.1
                    -π/4 # Negative 45 degrees
                elseif θ > π/2 || θ < -π/2
                    θ + π 
                elseif abs(θ - π/2) < 0.1 || abs(θ + π/2) < 0.1
                    π/2  # Ensure vertical text faces left
                else
                    θ
                end
            end

            if text
                
                if collinear

                    collinear_groups = collect(1:length(elements))
                    for i in 1:lastindex(elements)

                        typeof_i = typeof(elements[i]).parameters[1]

                        if typeof_i <: Union{Asap.FreeFree}
                            continue
                        end

                        for j in i+1:lastindex(elements)

                            typeof_j = typeof(elements[j]).parameters[1]

                            if typeof_j <: Union{Asap.FreeFree}
                                continue
                            end
                        
                            i_nodestart = elements[i].nodeStart.nodeID
                            i_nodeend = elements[i].nodeEnd.nodeID
                            j_nodestart = elements[j].nodeStart.nodeID
                            j_nodeend = elements[j].nodeEnd.nodeID

                            # Check each possible node pairing to find shared node
                            if i_nodestart == j_nodestart
                                if elements[i].nodeStart.dof[4:6] == [1,1,1] || elements[i].nodeStart.id == :column || typeof_i <: Union{Asap.FreeFixed, Asap.FreeFree} || typeof_j <: Union{Asap.FreeFixed, Asap.FreeFree}
                                    continue
                                end
                            elseif i_nodestart == j_nodeend
                                if elements[i].nodeStart.dof[4:6] == [1,1,1] || elements[i].nodeStart.id == :column || typeof_i <: Union{Asap.FreeFixed, Asap.FreeFree} || typeof_j <: Union{Asap.FixedFree, Asap.FreeFree}
                                    continue
                                end
                            elseif i_nodeend == j_nodestart
                                if elements[i].nodeEnd.dof[4:6] == [1,1,1] || elements[i].nodeEnd.id == :column || typeof_i <: Union{Asap.FixedFree, Asap.FreeFree} || typeof_j <: Union{Asap.FreeFixed, Asap.FreeFree}
                                    continue
                                end
                            elseif i_nodeend == j_nodeend
                                if elements[i].nodeEnd.dof[4:6] == [1,1,1] || elements[i].nodeEnd.id == :column || typeof_i <: Union{Asap.FixedFree, Asap.FreeFree} || typeof_j <: Union{Asap.FixedFree, Asap.FreeFree}
                                    continue
                                end
                            else
                                continue
                            end

                            if check_collinearity(elements[i], elements[j])
                                collinear_groups[j] = collinear_groups[i]
                            end
                        end
                    end

                    # Find the middle of the collinear group
                    group_indices = findall(x -> x == collinear_groups[i], collinear_groups)
                    mid_index = group_indices[Int(ceil(length(group_indices) / 2))]

                    group_length = sqrt((elements[group_indices[1]].nodeStart.position[1] - elements[group_indices[end]].nodeEnd.position[1])^2 + (elements[group_indices[1]].nodeStart.position[2] - elements[group_indices[end]].nodeEnd.position[2])^2)

                    if i == mid_index
                        
                        if iseven(length(group_indices))
                            # Determine font size based on beam length
                            base_fontsize = 8
                            multiplier = mini ? 1.5 : 2
                            fontsize = max(4, min(base_fontsize, element_length  * multiplier))

                            # Calculate midpoint coordinates for the middle beam
                            mid_x = (elements[mid_index].nodeStart.position[1] + elements[mid_index].nodeEnd.position[1]) / 2
                            mid_y = (elements[mid_index].nodeStart.position[2] + elements[mid_index].nodeEnd.position[2]) / 2
                        else
                            # Determine font size based on group length
                            base_fontsize = 8
                            multiplier = mini ? 1.5 : 2
                            fontsize = max(4, min(base_fontsize, group_length  * multiplier))

                            # Calculate midpoint coordinates for the collinear group
                            mid_x = mean([elements[k].nodeStart.position[1] + elements[k].nodeEnd.position[1] for k in group_indices]) / 2
                            mid_y = mean([elements[k].nodeStart.position[2] + elements[k].nodeEnd.position[2] for k in group_indices]) / 2
                        end

                        if !isempty(section_names)
                            push!(text_plots, (mid_x, mid_y, section_names[i], rotation, fontsize))
                        else
                            push!(text_plots, (mid_x, mid_y, sections[i], rotation, fontsize))
                        end

                    end
                else
                    # Determine font size based on beam length
                    base_fontsize = 8
                    multiplier = mini ? 1.5 : 2
                    fontsize = max(4, min(base_fontsize, element_length  * multiplier))

                    if !isempty(section_names)
                        push!(text_plots, (mid_x, mid_y, section_names[i], rotation, fontsize))
                    else
                        push!(text_plots, (mid_x, mid_y, sections[i], rotation, fontsize))
                    end
                end
            end

        end
        
    end

    for text_plot in text_plots
        mid_x, mid_y, section, rotation, fontsize = text_plot
        text!(new_ax, mid_x, mid_y, text=section,
            rotation=rotation,
            align=(:center, :center),
            fontsize=fontsize,
            color=:white,
            strokewidth=3,
            strokecolor=(:white, 0.8))
        text!(new_ax, mid_x, mid_y, text=section,
            rotation=rotation, 
            align=(:center, :center),
            fontsize=fontsize,
            color=:black)
    end

    for node in model.nodes[:column]
        scatter!(new_ax, node.position[1], node.position[2], color=:deeppink, markersize=10, marker=:rect)
    end

    GC.gc()

    return(new_fig)

end

function plot_slab(self::SlabAnalysisParams, results::Union{DataFrame,DataFrameRow}; text=true, mini=false, background=true, collinear=nothing)

    if !(results isa DataFrameRow)
        results = results[1, :] # Convert DataFrame to DataFrameRow by taking first row
    end

    sections = parse_sections(results.sections)

    if isnothing(collinear)
        collinear = results.collinear
    end

    return plot_slab(self, sections, text=text, mini=mini, background=background, collinear=collinear)

end

function plot_slab(self::SlabAnalysisParams, beam_sizing_params::SlabSizingParams; text=true, mini=false, background=true, collinear=nothing)
    
    slab_results = postprocess_slab(self, beam_sizing_params, check_collinear=false, concise=true);
    if occursin("W", slab_results.ids[1])
        sections = Vector{String}(slab_results.ids)
    else
        sections = parse.(Float64, slab_results.sections)
    end
    
    if isnothing(collinear)
        collinear = beam_sizing_params.collinear
    end

    return plot_slab(self, sections, text=text, mini=mini, background=background, collinear=collinear)

end

function plot_slab_section_delta(self::SlabAnalysisParams, delta_sections::Vector{Float64}; text=true, mini=false, background=true, collinear=false)

    model = self.model

    # Make a new figure and copy the axis
    if mini
        new_fig = Figure(size=(380,250))
    else
        new_fig = Figure(size=(600,400))
    end
    
    new_ax = Axis(new_fig[1,1], aspect=DataAspect())

    if background
        old_ax = self.plot_context.ax
        new_ax = copy_axis(old_ax, new_ax, alpha=0.5, scatter=false, beams=false)
    end

    hidespines!(new_ax)
    hidedecorations!(new_ax)
    
    # Copy relevant axis properties
    elements = self.model.elements[:beam]
    text_plots = []

    area_delta = delta_sections

    # Plot each beam element
    for (i, element) in enumerate(elements)
        if element.nodeStart.id == :wall && element.nodeEnd.id == :wall
            lines!(new_ax, [element.nodeStart.position[1], element.nodeEnd.position[1]], [element.nodeStart.position[2], element.nodeEnd.position[2]], color=:grey, linewidth=5)
            continue
        end
        
        x = [element.nodeStart.position[1], element.nodeEnd.position[1]]
        y = [element.nodeStart.position[2], element.nodeEnd.position[2]]
        linewidth = mini ? sqrt(abs(area_delta[i])) / 2 : sqrt(abs(area_delta[i]))
        color = area_delta[i]
        area_range = (minimum(area_delta), maximum(area_delta))

        # Calculate clipped line coordinates by moving inward from endpoints
        # Get release type from element and determine DOFs
        release_type = typeof(element).parameters[1]

        start_dof = if release_type <: Union{Asap.FixedFixed, Asap.FixedFree}  && sum(element.nodeStart.dof[4:6]) == 0
            [0,0,0,0,0,0] # Fixed start
        else
            [0,0,0,1,1,1] # Free start
        end
        end_dof = if release_type <: Union{Asap.FixedFixed, Asap.FreeFixed}  && sum(element.nodeEnd.dof[4:6]) == 0
            [0,0,0,0,0,0] # Fixed end
        else
            [0,0,0,1,1,1] # Free end
        end

        # Default to no clipping for moment connections
        if mini
            clip_amount = .2
        else
            clip_amount = .2
        end

        clip_start = sum(start_dof[4:6]) == 0 ? 0 : clip_amount  # Clip if rotation DOF is fixed
        clip_end = sum(end_dof[4:6]) == 0 ? 0 : clip_amount      # Clip if rotation DOF is fixed
        
        dx = x[2] - x[1]
        dy = y[2] - y[1]
        element_length = sqrt(dx^2 + dy^2)
        
        # Calculate parametric distances based on connection types
        t1 = clip_start/element_length
        t2 = (element_length-clip_end)/element_length
        
        x_clipped = [x[1] + t1*dx, x[1] + t2*dx]
        y_clipped = [y[1] + t1*dy, y[1] + t2*dy]
        lines!(new_ax, x_clipped, y_clipped, linewidth=linewidth, color=color, colorrange=area_range, colormap=:bam)

        # Calculate the length of the beam
        element_length = sqrt(dx^2 + dy^2)

        # Text
        # Calculate midpoint coordinates
        mid_x = (x[1] + x[2]) / 2
        mid_y = (y[1] + y[2]) / 2
        
        # Calculate rotation angle in radians
        rotation = let θ = atan((y[2] - y[1]), (x[2] - x[1]))
            if abs(θ - π/4) < 0.1
                π/4  # Positive 45 degrees
            elseif abs(θ + π/4) < 0.1
                -π/4 # Negative 45 degrees
            elseif θ > π/2 || θ < -π/2
                θ + π 
            elseif abs(θ - π/2) < 0.1 || abs(θ + π/2) < 0.1
                π/2  # Ensure vertical text faces left
            else
                θ
            end
        end

        if text
            
            if collinear

                collinear_groups = collect(1:length(elements))
                for i in 1:lastindex(elements)

                    typeof_i = typeof(elements[i]).parameters[1]

                    if typeof_i <: Union{Asap.FreeFree}
                        continue
                    end

                    for j in i+1:lastindex(elements)

                        typeof_j = typeof(elements[j]).parameters[1]

                        if typeof_j <: Union{Asap.FreeFree}
                            continue
                        end
                    
                        i_nodestart = elements[i].nodeStart.nodeID
                        i_nodeend = elements[i].nodeEnd.nodeID
                        j_nodestart = elements[j].nodeStart.nodeID
                        j_nodeend = elements[j].nodeEnd.nodeID

                        # Check each possible node pairing to find shared node
                        if i_nodestart == j_nodestart
                            if elements[i].nodeStart.dof[4:6] == [1,1,1] || elements[i].nodeStart.id == :column || typeof_i <: Union{Asap.FreeFixed, Asap.FreeFree} || typeof_j <: Union{Asap.FreeFixed, Asap.FreeFree}
                                continue
                            end
                        elseif i_nodestart == j_nodeend
                            if elements[i].nodeStart.dof[4:6] == [1,1,1] || elements[i].nodeStart.id == :column || typeof_i <: Union{Asap.FreeFixed, Asap.FreeFree} || typeof_j <: Union{Asap.FixedFree, Asap.FreeFree}
                                continue
                            end
                        elseif i_nodeend == j_nodestart
                            if elements[i].nodeEnd.dof[4:6] == [1,1,1] || elements[i].nodeEnd.id == :column || typeof_i <: Union{Asap.FixedFree, Asap.FreeFree} || typeof_j <: Union{Asap.FreeFixed, Asap.FreeFree}
                                continue
                            end
                        elseif i_nodeend == j_nodeend
                            if elements[i].nodeEnd.dof[4:6] == [1,1,1] || elements[i].nodeEnd.id == :column || typeof_i <: Union{Asap.FixedFree, Asap.FreeFree} || typeof_j <: Union{Asap.FixedFree, Asap.FreeFree}
                                continue
                            end
                        else
                            continue
                        end

                        if check_collinearity(elements[i], elements[j])
                            collinear_groups[j] = collinear_groups[i]
                        end
                    end
                end

                # Find the middle of the collinear group
                group_indices = findall(x -> x == collinear_groups[i], collinear_groups)
                mid_index = group_indices[Int(ceil(length(group_indices) / 2))]

                group_length = sqrt((elements[group_indices[1]].nodeStart.position[1] - elements[group_indices[end]].nodeEnd.position[1])^2 + (elements[group_indices[1]].nodeStart.position[2] - elements[group_indices[end]].nodeEnd.position[2])^2)

                if i == mid_index
                    
                    if iseven(length(group_indices))
                        # Determine font size based on beam length
                        base_fontsize = 8
                        multiplier = mini ? 1.5 : 2
                        fontsize = max(4, min(base_fontsize, element_length  * multiplier))

                        # Calculate midpoint coordinates for the middle beam
                        mid_x = (elements[mid_index].nodeStart.position[1] + elements[mid_index].nodeEnd.position[1]) / 2
                        mid_y = (elements[mid_index].nodeStart.position[2] + elements[mid_index].nodeEnd.position[2]) / 2
                    else
                        # Determine font size based on group length
                        base_fontsize = 8
                        multiplier = mini ? 1.5 : 2
                        fontsize = max(4, min(base_fontsize, group_length  * multiplier))

                        # Calculate midpoint coordinates for the collinear group
                        mid_x = mean([elements[k].nodeStart.position[1] + elements[k].nodeEnd.position[1] for k in group_indices]) / 2
                        mid_y = mean([elements[k].nodeStart.position[2] + elements[k].nodeEnd.position[2] for k in group_indices]) / 2
                    end

                    push!(text_plots, (mid_x, mid_y, area_delta[i], rotation, fontsize))
                end

            else
                # Determine font size based on beam length
                base_fontsize = 8
                multiplier = mini ? 1.5 : 2
                fontsize = max(4, min(base_fontsize, element_length  * multiplier))

                push!(text_plots, (mid_x, mid_y, area_delta[i], rotation, fontsize))
            end
        end
    end

    for text_plot in text_plots
        mid_x, mid_y, section, rotation, fontsize = text_plot
        text!(new_ax, mid_x, mid_y, text=string(round(section, digits=2)),
            rotation=rotation,
            align=(:center, :center),
            fontsize=fontsize,
            color=:white,
            strokewidth=3,
            strokecolor=(:white, 0.8))
        text!(new_ax, mid_x, mid_y, text=string(round(section, digits=2)),
            rotation=rotation, 
            align=(:center, :center),
            fontsize=fontsize,
            color=:black)
    end

    for node in model.nodes[:column]
        scatter!(new_ax, node.position[1], node.position[2], color=:deeppink, markersize=8, marker=:rect)
    end

    GC.gc()

    return(new_fig)

end
