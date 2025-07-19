function msfun_ica_meg_estimate(raw::Dict, data::Array, extdata=nothing, cfg=Dict())
    # Load required functions (assumed to be in other files or modules)
    using .MsFun  # assumes a module MsFun that includes these functions

    if extdata === nothing
        extdata = Array{Float64}(undef, 0)
    end

    shape = size(data)
    if length(shape) == 2
        epoching = false
        N, T = shape
        if !isempty(extdata) && (ndims(extdata) != 2 || size(extdata, 2) != T)
            error("Inconsistent extdata shape with data")
        end
        S = isempty(extdata) ? 0 : size(extdata, 1)
    elseif length(shape) == 3
        epoching = true
        K, N, T = shape
        if !isempty(extdata) && (ndims(extdata) != 3 || size(extdata, 1) != K || size(extdata, 3) != T)
            error("Inconsistent extdata shape with epoched data")
        end
        S = isempty(extdata) ? 0 : size(extdata, 2)
    else
        error("data must be 2D or 3D")
    end

    if N != 306
        println("msfun_ica_meg_estimate - WARNING: Data may not come from Neuromag Elekta MEG system")
    end

    println("msfun_ica_meg_estimate - Preparing data for ICA...")

    if epoching
        println("msfun_ica_meg_estimate - Baseline correcting and concatenating epochs...")
        data = data .- mean(data, dims=3)
        data = msfun_sig_concat_epoch(data, K, "epochnum")
    end

    IC = msfun_meg_ica_estimate(data, get(cfg, "ica", Dict()))
    IC = msfun_meg_ica_cumulantanalysis(IC, get(cfg, "cumulant", Dict()))

    if get(cfg, "corranalysis", false)
        println("msfun_ica_meg_estimate - Preparing external signals for correlation analysis...")
        if epoching
            println("msfun_ica_meg_estimate -       baseline correcting and concatenating epochs...")
            X = zeros(S, K*T)
            for k in 1:K
                av = mean(extdata[k, :, :], dims=3)
                for t in 1:T
                    X[:, (k-1)*T+t] = extdata[k, :, t] - vec(av)
                end
            end
            extdata = X
        end
        IC = msfun_meg_ica_corranalysis(IC, extdata, get(cfg, "corr", Dict()))
    end

    if epoching
        println("msfun_ica_meg_estimate - Restoring epochs in IC time courses...")
        IC["S"] = msfun_sig_concat_epoch(IC["S"], K, "epochlength")
    end

    if get(cfg, "spectralanalysis", false)
        IC = msfun_meg_ica_powerspectrum(IC, get(cfg, "fft", Dict()))
        IC = msfun_meg_ica_spectralanalysis(IC, get(cfg, "spectral", Dict()))
    end

    keep, reject, _ = msfun_meg_ica_resultfig(IC)
    if !haskey(cfg, "viewer")
        cfg["viewer"] = Dict()
    end
    cfg["viewer"]["list"] = reject

    msfun_meg_ica_viewer(IC, cfg["viewer"], raw)

    return IC
end
