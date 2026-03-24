include("_reinforced_concrete.jl")

b = 16.
d = 22.
l = ft_to_m(20.)
fc′ = 5.
fy = 60.

nodes = [Node([0.,0.,0.], [false, false, false, false, false, false]), Node([l,0.,0.], [true, false, false, true, true, true])]
section = Section(
    b * d * 0.0254^2,  # Area (A) in m^2 
    30e9,  # Elastic modulus (E) for concrete, in Pa
    12.5e9,   # Shear modulus (G) for concrete, in Pa
    in_to_m(b * d^3, factor=4) / 12,  # Second moment of area about x-axis (Ix) in m^4
    in_to_m(d * b^3, factor=4) / 12,  # Second moment of area about y-axis (Iy) in m^4
    in_to_m(b * d * (b^2 + d^2), factor=4) / 12  # Torsional constant (J) in m^4
)
elements = [Element(nodes[1], nodes[2], section)];
loads = [LineLoad(elements[1], [0,0,-kip_per_ft_to_N_per_m(9.4)])];
test_beam = Asap.Model(nodes, elements, loads);
Asap.solve!(test_beam);

internal_forces = InternalForces(elements[1], loads; resolution = 200);
beam_Mn = internal_forces.My .* (1/112.98) # Convert N⋅m to kip⋅in
beam_Vn = internal_forces.Vy .* (1/4448.22) # Convert N to kips
beam_x = in_to_m.(internal_forces.x,factor=-1)

maximum(abs.(beam_Mn))
maximum(abs.(beam_Vn))

lines(beam_Mn)

# Test select_rebar_gurobi function
zones, beam_ids = get_longitudinal_zones([beam_Mn], beam_x)
length(zones)

# Size the zones first
for zone in zones
    Mu = 0.9 * zone.max_Mn
    zone.As, zone.As′ = get_longitudinal_As_differentiable(b, d, 2.5, Mu, ReinforcedConcrete(5,60))
end
# Get tension steel
total_As_tension, per_zone_As_tension = select_rebar_gurobi(zones, beam_ids, b, clear_side=1.0, max_layers=2, max_bar_size=11, bars_per_layer=6, min_bar_number=2, max_unique_bars=2, force=:tension)
# Get compression steel 
total_As_compression, per_zone_As_compression = select_rebar_gurobi(zones, beam_ids, b, clear_side=1.0, max_layers=2, max_bar_size=11, bars_per_layer=6, min_bar_number=2, max_unique_bars=2, force=:compression)

print([zone.rebar_top for zone in zones])
print([zone.rebar_bottom for zone in zones])

println("Total tension steel area: ", total_As_tension)
println("Total compression steel area: ", total_As_compression)



required_As, required_As′, required_Av, h, d′ = size_concrete_section(b, d, beam_x, beam_Mn, beam_Vn, fc′, fy, d′=2.5, rounded=false, α=π/2, plot=true);
required_As, required_As′, required_Av, h, d′ = size_concrete_section(10., 17.5, beam_x, [1600/0.9], beam_Vn, 4., fy, d′=2.5, rounded=false, plot=true);

r = size_concrete_zones(M_zones) 
println(r.minimizer)
println(r.minimum)

lb = [10.0, 10.0]
ub = [25.0, 30.0]
x0 = [15.0, 15.0]
test = [15.0, 20.5]

objective_function, constraint_function = get_rc_objective(M_zones, ReinforcedConcrete(5,60), max_h=24.0, plot=false);
objective_function(test)
constraint_function(test)

println(required_As, ", ", required_As′, ", ", required_Av, ", ", h, ", ", d′)

model = Nonconvex.Model(objective_function);
Nonconvex.addvar!(model, lb, ub);
Nonconvex.add_ineq_constraint!(model, x -> constraint_function(x));
# aspect ratio
Nonconvex.add_ineq_constraint!(model, x -> x[2]/x[1] - 2.5); # d/b <= 2.5
Nonconvex.add_ineq_constraint!(model, x -> 1.0 - x[2]/x[1]); # d/b >= 1.0

alg = NLoptAlg(:LN_COBYLA)
options = NLoptOptions()
r = optimize(model, alg, x0, options = options)
r.minimizer
r.minimum

required_As, required_As′, required_Av, h, d′ = size_concrete_section(r.minimizer[1], r.minimizer[2], beam_x, beam_Mn, beam_Vn, fc′, fy, d′=2.5, rounded=false, α=π/4, plot=true);

"""
FINITE ELEMENT METHOD
- section modulus of bh2/6
- assume T = C strength
- 2% reinforcement ratio

TODO
- implement zonal approach
- implement development length
"""

