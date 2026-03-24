using AsapOptim, Asap, AsapToolkit
using Zygote
using LinearSolve, LinearAlgebra

begin
    w_section = W("W460X158")
    @show w_section.name

    section = Section(
        Steel_kNm,
        w_section.A / 1e6,
        w_section.Ix / 1e12,
        w_section.Iy / 1e12,
        w_section.J / 1e12
    )
end;

# generate
n = 14
begin
    nx = n
    dx = 1.5
    ny = n
    dy = 1.75
    dz = 2.

    # loads
    load = [0., 0., -20]

    spaceframe = SpaceFrame(nx, dx, ny, dy, dz, section; load = load)
    model = spaceframe.model
end

begin
    @assert n % 2 == 0

    igrid = spaceframe.itop

    imid = Int(n / 2)

    iparent = igrid[2:imid, 2:imid]

    ichild1 = reverse(igrid[2:imid, imid+1:end-1], dims = 2)
    factors1 = [-1., 1.]

    ichild2 = reverse(igrid[imid+1:end-1, 2:imid], dims = 1)
    factors2 = [1., -1.]

    ichild3 = reverse(igrid[imid+1:end-1, imid+1:end-1])
    factors3 = [-1., -1.]

    fac = 0.9
    x = dx * fac / 2
    y = dy * fac / 2
    z = 1.5
end

begin
    # make variables
    vars = Vector{TrussVariable}()
    # coupled_vars = FrameVariable[]

    fac = .9
    x = dx * fac / 2
    y = dy * fac / 2
    z = 1.5


    for i in eachindex(iparent)

        i0 = iparent[i]
        i1 = ichild1[i]
        i2 = ichild2[i]
        i3 = ichild3[i]

        push!(vars, SpatialVariable(i0, 0., -x, x, :X))
        push!(vars, SpatialVariable(i1, 0., -x, x, :X))
        push!(vars, SpatialVariable(i2, 0., -x, x, :X))
        push!(vars, SpatialVariable(i3, 0., -x, x, :X))

        push!(vars, SpatialVariable(i0, 0., -y, y, :Y))
        push!(vars, SpatialVariable(i1, 0., -y, y, :Y))
        push!(vars, SpatialVariable(i2, 0., -y, y, :Y))
        push!(vars, SpatialVariable(i3, 0., -y, y, :Y))
    end

    # for i in igrid
    #     push!(vars, SpatialVariable(i, 0., -x, x, :X))
    #     push!(vars, SpatialVariable(i, 0., -y, y, :Y))
    # end
end

# vars = TrussVariable[
#     [SpatialVariable(i, 0., -x, x, :X) for i in igrid][:];
#     [SpatialVariable(i, 0., -y, y, :Y) for i in igrid][:]
# ]

params = TrussOptParams(model, vars)