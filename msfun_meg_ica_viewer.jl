# msfun_meg_ica_viewer.jl

"""
msfun_meg_ica_viewer(raw, times, IC; cfg=Dict())

This is a placeholder function for the MATLAB `meg_ica_viewer`.
It was originally an interactive GUI viewer for inspecting IC time series and topoplots.

Arguments:
- `raw`: header information or data structure
- `times`: time vector or matrix
- `IC`: dictionary with fields "S", "A", possibly "powspctrm"
- `cfg`: configuration dictionary (optional)

Note: Interactive GUI not implemented in this Julia version.
"""

function msfun_meg_ica_viewer(raw, times, IC; cfg=Dict())
    error("Interactive GUI viewer must be rebuilt using Makie.jl, Gtk.jl, or similar in Julia.")
end