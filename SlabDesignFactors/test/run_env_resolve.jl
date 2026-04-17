#!/usr/bin/env julia
# Environment sync for the *currently running* Julia (e.g. cluster `module load julia/1.10.4`).
#
# Same idea as `SlabDesignFactors/scripts/_scripts.jl` + `run_preflight.slurm`: activate the
# repo project and `Pkg.instantiate()` using **Manifest.toml**. Do not delete Manifest and call
# bare `Pkg.resolve()` — AsapOptim / AsapToolkit are **path** packages, not General-registry
# names; resolve-from-Project.toml-only would error with "expected package ... to be registered".
#
# Usage (from repo root):
#   julia --project=. SlabDesignFactors/test/run_env_resolve.jl
#
# Optional — force rebuilding Manifest.toml.
# Also, if Manifest.toml is missing, this script automatically bootstraps path deps and
# resolves a fresh manifest even when REBUILD_MANIFEST is not set.
# Julia 1.10’s Pkg does not treat Project.toml [sources] the same as newer versions for
# bare `resolve()`.
#   REBUILD_MANIFEST=1 julia --project=. SlabDesignFactors/test/run_env_resolve.jl

using Pkg

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PROJ = joinpath(ROOT, "Project.toml")
const MANIFEST = joinpath(ROOT, "Manifest.toml")
const MANIFEST_BAK = joinpath(ROOT, "Manifest.toml.before_rebuild")

isfile(PROJ) || error("Expected Project.toml at $PROJ (wrong cwd?)")

println("Julia VERSION: ", VERSION)
println("Project root: ", ROOT)
println("REBUILD_MANIFEST: ", repr(get(ENV, "REBUILD_MANIFEST", "")))

Pkg.activate(ROOT)
println("Active project: ", Base.active_project())

function _develop_local_path_deps!()
    # Path packages: must be developed before `resolve` on Julia 1.10
    # (General registry has no AsapToolkit / AsapOptim).
    asap_toolkit = joinpath(ROOT, "AsapToolkit")
    asap_optim = joinpath(ROOT, "AsapOptim")
    isdir(asap_toolkit) || error("Missing path package directory: $asap_toolkit")
    isdir(asap_optim) || error("Missing path package directory: $asap_optim")
    println("==> Pkg.develop — AsapToolkit, AsapOptim (local paths)")
    Pkg.develop(Pkg.PackageSpec(path=asap_toolkit))
    Pkg.develop(Pkg.PackageSpec(path=asap_optim))
end

needs_rebuild = (get(ENV, "REBUILD_MANIFEST", "") == "1") || !isfile(MANIFEST)

if needs_rebuild
    if isfile(MANIFEST)
        println("==> REBUILD_MANIFEST=1: backing up Manifest.toml → $(basename(MANIFEST_BAK))")
        cp(MANIFEST, MANIFEST_BAK; force=true)
        rm(MANIFEST)
    else
        println("==> Manifest.toml missing: bootstrapping a fresh manifest")
    end
    _develop_local_path_deps!()
    println("==> Pkg.resolve()")
    Pkg.resolve()
    println("    resolve: OK")
end

println("==> Pkg.instantiate()")
Pkg.instantiate()
println("    instantiate: OK")

println("==> Pkg.precompile()")
Pkg.precompile()
println("    precompile: OK")

println("==> Smoke: core table stack (Tables + DataFrames + CSV)")
using Tables, DataFrames, CSV
println("    smoke: OK")

println("==> Pkg.status()")
Pkg.status()

println("--- run_env_resolve.jl finished ---")
