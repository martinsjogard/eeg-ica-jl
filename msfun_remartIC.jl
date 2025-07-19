using LinearAlgebra
using Plots

function msfun_remartIC(data, IC, lst=nothing, flag=true)
    if lst === nothing || isempty(lst)
        if haskey(IC, :artdetect) && haskey(IC[:artdetect], :list)
            lst = IC[:artdetect][:list]
        else
            error("List of artifactual ICs must be provided or found in IC[:artdetect][:list].")
        end
    end

    data2 = copy(data)

    if ndims(data) == 2
        A, S = IC[:A], IC[:S]
        if size(data, 1) != size(A, 1) || size(data, 2) != size(S, 2) || size(S, 1) != size(A, 2)
            error("Input data and IC structure inconsistent for continuous data.")
        end
        data2 = data - A[:, lst] * S[lst, :]
    elseif ndims(data) == 3
        A, S = IC[:A], IC[:S]
        if size(data, 2) != size(A, 1) || size(data, 1) != size(S, 1) || size(data, 3) != size(S, 3) || size(S, 2) != size(A, 2)
            error("Input data and IC structure inconsistent for epoched data.")
        end
        for k in 1:size(data, 1)
            data2[k, :, :] = data[k, :, :] .- A[:, lst] * S[k, lst, :]
        end
    else
        error("Data must be a 2D or 3D array.")
    end

    if flag
        nfig = length(lst)
        plot(layout=(nfig, 1), size=(800, 200 * nfig))
        if ndims(data) == 2
            for (i, ic) in enumerate(lst)
                chan = argmax(abs.(IC[:A][:, ic]))
                p = plot(data[chan, :], label="Original")
                plot!(p, data2[chan, :], label="Cleaned", color=:red)
                title!(p, "Channel $chan for IC $ic")
                plot!(p, legend=false, framestyle=:none)
                plot!(p)
            end
        else
            for (i, ic) in enumerate(lst)
                chan = argmax(abs.(IC[:A][:, ic]))
                diffs = sum((data[:, chan, :] .- IC[:S][:, ic, :]).^2, dims=2)
                trial = argmax(vec(diffs))
                p = plot(vec(data[trial, chan, :]), label="Original")
                plot!(p, vec(data2[trial, chan, :]), label="Cleaned", color=:red)
                title!(p, "Channel $chan, Trial $trial for IC $ic")
                plot!(p, legend=false, framestyle=:none)
                plot!(p)
            end
        end
    end

    return data2
end