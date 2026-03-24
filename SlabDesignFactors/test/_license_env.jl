# Use repo-local Gurobi license before Gurobi is loaded (mirrors `env_gurobi.sh`).
# Safe to `include` multiple times.
let repo = normpath(joinpath(@__DIR__, "..", "..")),
    lic = joinpath(repo, "secrets", "gurobi.lic")

    if isfile(lic)
        ENV["GRB_LICENSE_FILE"] = lic
    end
end
