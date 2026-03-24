# load dependencies
using Pkg; Pkg.activate(".")
using AsapOptim
import AsapOptim: Asap

#=
CURRENT MWE FOR SEGFAULT

Repeatedly execute everything below starting at `n = ...` in the REPL (i.e. press shift + enter).
If you reach the end of the script, try changing the value of n to some even integer value (e.g. n = 26, 40, 30, 18, ...) and rerun the script.
Eventually, a fatal crash (segfault) will occur. This may happen on the first time you run this script, or the fifth, etc.

Notes:
- Generally, this occurs sooner with larger values of n.
- Occurs in all version of Julia 1.9 and higher (this environment was made in 1.10.3). DOES NOT OCCUR in the latest version of Julia 1.8 (1.8.5).
- Occurs consistently on a M1 Pro Macbook Pro and a Windows 11 Intel PC.
- This does NOT occur if you enclose everything in a for loop and run for a bunch of iterations
=#


# select some even integer value for n (number of nodes in the square grid)
# hyperparameters
begin
    # structural cross section assigned to all elements
    section = Asap.Section(
        20e-3, # A [m²]
        2e8, # E [kN/m²]
        8e7, # G [kN/m²]
        7.95e-4, # Ix [m⁴]
        9.2e-5, # Iy [m⁴]
        3.11e-6 # J [m⁴]
    )

    # grid dimensions
    Lx = 20 #m
    Ly = 20 #m
end

#=
=#

n = 16
begin
    # nodal positions
    x_positions = range(0, Lx, n)
    y_positions = range(0, Ly, n)

    # distance between nodes
    dx = Lx / (n-1)
    dy = Ly / (n-1)

    # node positions
    xyz = Vector{Vector{Float64}}()
    Xmatrix = zeros(Float64, n, n)
    Ymatrix = zeros(Float64, n, n)
    igrid = zeros(Int64, n, n)

    index = 1
    for iy = 1:n
        for ix = 1:n

            x = x_positions[iy]
            y = y_positions[ix]

            igrid[ix, iy] = index
            index += 1

            push!(xyz, [x, y, 0.])
            Xmatrix[ix, iy] = x
            Ymatrix[ix, iy] = y

        end
    end

    # extract corner node indices to fix
    support_indices = [igrid[1, 1], igrid[n, 1], igrid[1, n], igrid[n, n]]

    # [Asap] fix the corner indices as supports
    nodes = [Asap.Node(pos, :free, :free) for pos in xyz]

    # [Asap] make support nodes
    for node in nodes[support_indices]
        Asap.fixnode!(node, :pinned)
        node.id = :support
    end

    # [Asap] element collector
    elements = Vector{Asap.Element}()

    # [Asap] horizontal elements
    for i = 1:n
        for j = 1:n-1
            index = [igrid[i,j], igrid[i,j+1]]
            push!(elements, Asap.Element(nodes[index]..., section))
        end
    end

    # [Asap] vertical elements
    for j = 1:n
        for i = 1:n-1
            index = [igrid[i,j], igrid[i+1,j]]
            push!(elements, Asap.Element(nodes[index]..., section)) 
        end
    end

    # [Asap] loads
    # applied load vector
    load = [0., 0., -20] #kN
    loads = [Asap.NodeForce(node, load) for node in nodes[:free]]

    # [Asap] assemble model and solve
    model = Asap.Model(nodes, elements, loads)
    Asap.solve!(model)
end;

begin
    igrid = reshape(1:n^2, n, n)

    imid = Int(floor(n / 2))

    iparent = igrid[2:imid, 2:imid]

    ichild1 = reverse(igrid[2:imid, imid+1:2imid-1], dims = 2)
    factors1 = [-1., 1.]

    ichild2 = reverse(igrid[imid+1:2imid-1, 2:imid], dims = 1)
    factors2 = [1., -1.]

    ichild3 = reverse(igrid[imid+1:2imid-1, imid+1:2imid-1])
    factors3 = [-1., -1.]
end

# [AsapOptim] make design variables
begin
    vars = Vector{FrameVariable}()

    for i in eachindex(iparent)

        i0 = iparent[i]
        i1 = ichild1[i]
        i2 = ichild2[i]
        i3 = ichild3[i]

        push!(vars, SpatialVariable(i0, 0., -dx, dx, :X))
        push!(vars, SpatialVariable(i1, 0., -dx, dx, :X))
        push!(vars, SpatialVariable(i2, 0., -dx, dx, :X))
        push!(vars, SpatialVariable(i3, 0., -dx, dx, :X))

        push!(vars, SpatialVariable(i0, 0., -dy, dy, :Y))
        push!(vars, SpatialVariable(i1, 0., -dy, dy, :Y))
        push!(vars, SpatialVariable(i2, 0., -dy, dy, :Y))
        push!(vars, SpatialVariable(i3, 0., -dy, dy, :Y))
    end

    AsapOptim.process_variables!(vars)
end

# begin
#     # make variables
#     vars = Vector{FrameVariable}()

#     for i in igrid[2:end-1, 2:end-1]
#         push!(vars, SpatialVariable(i, 0., -dx, dx, :X))
#         push!(vars, SpatialVariable(i, 0., -dy, dy, :Y))
#         push!(vars, SpatialVariable(i, 0., -dy, dy, :Z))
#     end


#     vals, l, u  = AsapOptim.process_variables!(vars)

# end

# [AsapOptim] make optimization parameters
# FrameOptParams2(model, vars)

#=
If you've come this far without a segfault, congrats! Trying changing the value of `n` above or changing some hyperparameters (Lx, x, y, ....) and rerunning.
=#