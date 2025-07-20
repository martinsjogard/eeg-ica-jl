using LinearAlgebra
using Plots
using Statistics

function msfun_remartIC(data, IC::Dict, list=nothing, flag=true)
    if list === nothing || isempty(list)
        if haskey(IC, "artdetect") && haskey(IC["artdetect"], "list")
            list = IC["artdetect"]["list"]
        else
            error("remartIC - ERROR : A list of artifactual ICs must be given either in IC.artdetect.list or in a third input vector list...")
        end
    end

    data2 = copy(data)
    println("remartIC - Removing $(length(list)) ICs...")

    if ndims(data) == 2
        data2 .= data .- IC["A"][:, list] * IC["S"][list, :]
    elseif ndims(data) == 3
        for k in 1:size(data, 1)
            data2[k, :, :] .= data[k, :, :] .- IC["A"][:, list] * IC["S"][k, list, :]
        end
    else
        error("remartIC - ERROR : First input data must be array with 2 or 3 dimensions...")
    end

    if flag
        println("remartIC - Generating comparative plots...")
        nfig = length(list)
        plot_list = []

        if ndims(data) == 2
            for (idx, ic) in enumerate(list)
                chan = argmax(abs.(IC["A"][:, ic]))
                p = plot(data[chan, :], label="Original")
                plot!(p, data2[chan, :], label="Cleaned", color=:red)
                title!(p, "Channel $(chan) with max mixing for IC $(ic)")
                push!(plot_list, p)
            end
        else
            for (idx, ic) in enumerate(list)
                chan = argmax(abs.(IC["A"][:, ic]))
                diffs = mapslices(x -> sum(abs2, x), data[:, chan, :] .- IC["S"][:, ic, :], dims=2)
                trial = argmax(diffs)
                p = plot(vec(data[trial, chan, :]), label="Original")
                plot!(p, vec(data2[trial, chan, :]), label="Cleaned", color=:red)
                title!(p, "Channel $(chan), trial $(trial) with max mixing for IC $(ic)")
                push!(plot_list, p)
            end
        end
        display(plot(plot_list..., layout=(nfig, 1)))
    end

    println("remartIC - ICs removed from data.")
    return data2
end