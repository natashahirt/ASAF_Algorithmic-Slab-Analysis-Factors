"""
Justifying As = Mu/4d 
source spreadsheet: https://docs.google.com/spreadsheets/d/0B2YoJJbAdgUzVUtNY1FMZHJNMHc/edit?
"""

function ρ_estimate_Mu4d(fc′, fy)
    b = d = 1    
    β1 = get_β1(fc′)
    As_min = singly_reinforced_min_ρ(fc′, fy) * (b * d)
    As_max = singly_reinforced_max_ρ(fc′, fy) * (b * d)
    As_avg = (As_min + As_max) / 2
    z_min = 1-0.5*(As_max*fy)/(β1*fc′*1) # min lever arm ratio
    z_max = 1-0.5*(As_min*fy)/(β1*fc′*1) # max lever arm ratio
    z_avg = (z_min + z_max) / 2
    ϕMn = As_max * z_avg * fy / 12 * 0.9
    coefficient = ϕMn / As_max
    return coefficient
end

ρ_estimate_Mu4d(4,60)

begin
    fig = Figure();
    ax = Axis(fig[1,1], 
        xlabel = "Steel Strength fy (ksi)", 
        ylabel = "Average ϕMn/As (kip-in/in²)",
        title = "Moment Capacity per Unit Steel Area")

    # Use the exact x points from fy_vec
    fy_points = [0, 5, 10, 30, 40, 50, 60, 70, 80, 120, 160, 180, 200]
    fc_vec = [3,4,5,6,7,8]
    for (i,fy) in enumerate(fy_points)
        coefficient_mean = minimum([ρ_estimate_Mu4d(fc′, fy) for fc′ in fc_vec])
        scatter!(ax, i, coefficient_mean)
    end

    ax.xticks = (1:length(fy_points), string.(fy_points))

    display(fig)
end

begin
    fig = Figure();
    ax = Axis(fig[1,1], 
        xlabel = "Steel Strength fy (ksi)", 
        ylabel = "Average ϕMn/As (kip-in/in²)",
        title = "Moment Capacity per Unit Steel Area")

    # Use the exact x points from fy_vec
    fy_points = range(0, stop=200, step=10)
    fc_vec = [3,4,5,6,7,8]
    for (i,fy) in enumerate(fy_points)
        coefficient_mean = minimum([ρ_estimate_Mu4d(fc′, fy) for fc′ in fc_vec])
        scatter!(ax, i, coefficient_mean)
    end

    display(fig)
end