using Statistics, Printf, LinearAlgebra, Plots

function msfun_ica_meg_plot(IC::Dict)
    reject = Int[]

    Ncumulant = 0
    if haskey(IC, "cumulant")
        Ncumulant += haskey(IC["cumulant"], "skew") + haskey(IC["cumulant"], "kurt")
        append!(reject, get(IC["cumulant"], "list", Int[]))
    end

    corrnames = String[]
    Ncorr = 0
    if haskey(IC, "corr")
        corrnames = filter(x -> !(x in ["list", "Tcorr"]), keys(IC["corr"]))
        Ncorr = length(corrnames)
        append!(reject, get(IC["corr"], "list", Int[]))
    end

    Nspectral = 0
    if haskey(IC, "spectral") && haskey(IC["spectral"], "gof")
        Nspectral = size(IC["spectral"]["gof"], 2)
        append!(reject, get(IC["spectral"], "list", Int[]))
    end

    reject = sort(unique(reject))
    numofic = size(IC["A"], 2)
    keep = setdiff(1:numofic, reject)

    Ntot = Ncumulant + Ncorr + Nspectral
    ncol = floor(Int, sqrt(Ntot))
    nrow = ceil(Int, Ntot / ncol)

    println("msfun_ica_meg_plot - Generating plots...")
    plot_list = []

    # Plot skewness
    if haskey(IC, "cumulant")
        if haskey(IC["cumulant"], "skew")
            p = plot(IC["cumulant"]["skew"], seriestype=:scatter, label="", title="CUMULANT ANALYSIS : SKEWNESS")
            scatter!(p, reject, IC["cumulant"]["skew"][reject], label="", color=:red)
            hline!(p, [IC["cumulant"]["Tskew"], -IC["cumulant"]["Tskew"]], linestyle=:dash, color=:green, label="")
            xlabel!(p, "IC index"); ylabel!(p, "skew")
            push!(plot_list, p)
        end
        if haskey(IC["cumulant"], "kurt")
            p = plot(IC["cumulant"]["kurt"], seriestype=:scatter, label="", title="CUMULANT ANALYSIS : KURTOSIS")
            scatter!(p, reject, IC["cumulant"]["kurt"][reject], label="", color=:red)
            hline!(p, [IC["cumulant"]["Tkurt"]], linestyle=:dash, color=:green, label="")
            xlabel!(p, "IC index"); ylabel!(p, "kurt")
            push!(plot_list, p)
        end
    end

    for name in corrnames
        v = IC["corr"][name]
        Tcorr = IC["corr"]["Tcorr"]
        p = plot(v, seriestype=:scatter, label="", title="CORR ANALYSIS : IC/$(name)")
        scatter!(p, reject, v[reject], label="", color=:red)
        hline!(p, [Tcorr, -Tcorr], linestyle=:dash, color=:green, label="")
        xlabel!(p, "IC index"); ylabel!(p, "corr")
        push!(plot_list, p)
    end

    for k in 1:Nspectral
        gof_vals = IC["spectral"]["gof"][:, k]
        threshold = IC["spectral"]["Tgof"][k]
        label = IC["spectral"]["fit"][2k - 1]
        range_vals = IC["spectral"]["fit"][2k]
        p = plot(gof_vals, seriestype=:scatter, label="",
                 title="SPECTRAL ANALYSIS : $(label) on [$(range_vals[1]),$(range_vals[2])]")
        scatter!(p, reject, gof_vals[reject], label="", color=:red)
        hline!(p, [threshold], linestyle=:dash, color=:green, label="")
        xlabel!(p, "IC index"); ylabel!(p, "gof")
        push!(plot_list, p)
    end

    fig = plot(plot_list..., layout=(nrow, ncol), size=(400*ncol, 300*nrow))
    return keep, reject, fig
end
