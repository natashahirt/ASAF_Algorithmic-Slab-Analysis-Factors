"""
Headless / no-Makie mode for cluster nodes (Linux without DISPLAY/WAYLAND) or when
system `glib` conflicts with JLLs.

Set `SLABDESIGN_SKIP_MAKIE=1` to force skip, or `=0` to force CairoMakie (e.g. X11 forwarding).
"""
function _slab_skip_makie()
    if haskey(ENV, "SLABDESIGN_SKIP_MAKIE")
        v = strip(get(ENV, "SLABDESIGN_SKIP_MAKIE", ""))
        isempty(v) && return _slab_infer_headless_skip_makie()
        lv = lowercase(v)
        if lv in ("0", "false", "no", "n", "off")
            return false
        end
        return lv in ("1", "true", "yes", "y", "on")
    end
    return _slab_infer_headless_skip_makie()
end

function _slab_infer_headless_skip_makie()
    Sys.iswindows() && return false
    Sys.isapple() && return false
    return isempty(get(ENV, "DISPLAY", "")) && isempty(get(ENV, "WAYLAND_DISPLAY", ""))
end

"""Return whether to skip Makie; sync `ENV[\"SLABDESIGN_SKIP_MAKIE\"]` to `\"0\"` or `\"1\"`."""
function _slab_sync_skip_makie_env!()
    skip = _slab_skip_makie()
    ENV["SLABDESIGN_SKIP_MAKIE"] = skip ? "1" : "0"
    return skip
end
