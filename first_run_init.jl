using Pkg

# Regenerate Manifest.toml from Project.toml (run with the Julia version you want recorded).
# Example:
#   /path/to/julia/1.10.4/bin/julia --project=. first_run_init.jl
#
# Path packages must be linked before instantiate/resolve.
#
# Precompile: on Windows and macOS, full Pkg.precompile() is the default. On Linux with no
# DISPLAY/WAYLAND, we assume a headless cluster login node and skip Makie/CairoMakie (override
# with SKIP_VIZ_PRECOMPILE=0 to force full precompile, or =1 to force skip anywhere).

const ROOT = dirname(@__FILE__)

"""
Whether to skip Makie/Cairo precompile. If `SKIP_VIZ_PRECOMPILE` is set to a non-empty value,
it wins (0/false/no = full precompile). Otherwise we infer: Windows/macOS → full precompile;
Linux/FreeBSD with no DISPLAY and no WAYLAND_DISPLAY → skip (typical SSH / batch node).
"""
function _skip_viz_precompile()
    if haskey(ENV, "SKIP_VIZ_PRECOMPILE")
        v = strip(get(ENV, "SKIP_VIZ_PRECOMPILE", ""))
        if isempty(v)
            return _infer_headless_skip_viz()
        end
        lv = lowercase(v)
        if lv in ("0", "false", "no", "n", "off")
            return false
        end
        return lv in ("1", "true", "yes", "y", "on")
    end
    return _infer_headless_skip_viz()
end

function _infer_headless_skip_viz()
    Sys.iswindows() && return false
    Sys.isapple() && return false
    # Linux / FreeBSD: no display server → typical HPC login / job node without GUI.
    return isempty(get(ENV, "DISPLAY", "")) && isempty(get(ENV, "WAYLAND_DISPLAY", ""))
end

"""Packages to omit from precompile when skipping the Makie/Cairo stack."""
function _skip_viz_pkg(name::AbstractString)
    name == "Cairo" && return true
    return occursin("Makie", String(name))
end

"""Full `Pkg.precompile()`, or all deps except the Makie/Cairo stack when `skip_viz` is true."""
function _precompile_for_environment(skip_viz::Bool)
    if !skip_viz
        Pkg.precompile()
        return
    end
    @info "Skipping Makie/CairoMakie precompile (headless or SKIP_VIZ_PRECOMPILE): precompiling remaining packages only."
    names = String[]
    for (_, p) in Pkg.dependencies()
        p.name === nothing && continue
        _skip_viz_pkg(p.name) && continue
        push!(names, p.name)
    end
    sort!(unique!(names))
    isempty(names) && return
    Pkg.precompile(names)
end

skip_viz = _skip_viz_precompile()

Pkg.activate(ROOT)

# Asap resolves from the General registry (Project.toml).

Pkg.develop(path=joinpath(ROOT, "AsapToolkit"))
Pkg.develop(path=joinpath(ROOT, "AsapOptim"))
if skip_viz
    # Otherwise `instantiate` may precompile CairoMakie/Makie before we run selective precompile.
    ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
    @info "JULIA_PKG_PRECOMPILE_AUTO=0 for instantiate (avoids precompiling viz stack during Pkg.instantiate)."
end
Pkg.instantiate()
# Optional:
# Pkg.build("Gurobi")

# Fewer parallel precompile workers avoids Makie/Cairo pidfile races on login nodes.
ENV["JULIA_NUM_PRECOMPILE_TASKS"] = get(ENV, "JULIA_NUM_PRECOMPILE_TASKS", "1")
_precompile_for_environment(skip_viz)