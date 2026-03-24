using Asap, AsapOptim

# frame Optimization
begin
    section = Section(
        20e-3,
        2e8,
        8e7,
        7.95e-4,
        9.2e-5,
        3.11e-6
    )
end

begin
    Lx = 25.
    Ly = 15.
    x = 1.25
    y = 1.45
    z = 2.5
end

n = 18

# generate
begin
    support = :corner
    support_type = :pinned

    # loads
    load = [0., 0., -20]

    x_positions = range(0, Lx, n)
    y_positions = range(0, Ly, n)

    dx = Lx / (n-1)
    dy = Ly / (n-1)

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

    if support == :corner
        support_indices = [igrid[1, 1], igrid[n, 1], igrid[1, n], igrid[n, n]]
    elseif support == :x
        support_indices = igrid[[1, n], :]
    elseif support == :y
        support_indices = igrid[:, [1, n]]
    else
        support_indices = [igrid[[1, n], :][:]; igrid[2:n-1, [1, n]][:]]
    end

    #make nodes
    nodes = [Node(pos, :free, :free) for pos in xyz]

    #make support nodes
    for node in nodes[support_indices]
        fixnode!(node, support_type)
        node.id = :support
    end

    #make elements
    elements = Vector{Element}()

    #horizontal elements
    for i = 1:n
        for j = 1:n-1
            index = [igrid[i,j], igrid[i,j+1]]
            push!(elements, Element(nodes[index]..., section))
        end
    end

    #vertical elements
    for j = 1:n
        for i = 1:n-1
            index = [igrid[i,j], igrid[i+1,j]]
            push!(elements, Element(nodes[index]..., section)) 
        end
    end

    #loads
    loads = [NodeForce(node, load) for node in nodes[:free]]

    #assemble
    model = Model(nodes, elements, loads)
    Asap.solve!(model)
end;

# design variables
begin
    # n = 30
    @assert n % 2 == 0

    igrid = reshape(1:n^2, n, n)

    imid = Int(n / 2)

    iparent = igrid[2:imid, 2:imid]

    ichild1 = reverse(igrid[2:imid, imid+1:end-1], dims = 2)
    factors1 = [-1., 1.]

    ichild2 = reverse(igrid[imid+1:end-1, 2:imid], dims = 1)
    factors2 = [1., -1.]

    ichild3 = reverse(igrid[imid+1:end-1, imid+1:end-1])
    factors3 = [-1., -1.]
end

begin
    # make variables
    vars = Vector{FrameVariable}()
    # vars = []

    # fac = .9
    # x = dx * fac / 2
    # y = dy * fac / 2
    # z = 1.5


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
end

params = FrameOptParams(model, vars)

x0 = deepcopy(params.values)

@time res = solve_frame(x0, params);