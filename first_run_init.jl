using Pkg

# Activate the current environment
Pkg.activate(".")

# Add necessary packages
Pkg.add("Asap")
Pkg.add("AsapOptim")
Pkg.add(url="https://github.com/keithjlee/AsapToolkit")
Pkg.add("CairoMakie")
Pkg.add("ChainRulesCore")
Pkg.add("Colors")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("GLMakie")
Pkg.add("HTTP")
Pkg.add("Interpolations")
Pkg.add("JSON")
Pkg.add("LinearSolve");
Pkg.add("Nonconvex")
Pkg.add("Revise")
Pkg.add("Sockets")
Pkg.add("Statistics")
Pkg.add("StatsBase")
Pkg.add("UnPack")
Pkg.add("Zygote")

# Redirect to local directories
Pkg.develop(path="./AsapToolkit")
Pkg.develop(path="./AsapOptim")

# Resolve and precompile packages
Pkg.resolve()
Pkg.precompile()

# Instantiate the project
Pkg.instantiate()