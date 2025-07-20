using FFTW
using Statistics
using LinearAlgebra

include("msfun_filt_preparecosine.jl")

function msfun_ica_meg_signalcorrestimate_trialwise(IC::Dict, extdata::Array, cfg::Dict)
    if extdata === nothing || isempty(extdata)
        println("msfun_ica_meg_signalcorrestimate_trialwise - WARNING : No external data supplied... Skipping the correlation analysis.")
        return IC
    end

    if !haskey(IC, "S") || !(IC["S"] isa AbstractArray) || !(ndims(IC["S"]) in (2, 3))
        error("msfun_ica_meg_signalcorrestimate_trialwise - ERROR : IC structure missing elements or inconsistent.")
    end

    epoching = ndims(IC["S"]) == 3
    if epoching
        K, numofic, T = size(IC["S"])
    else
        numofic, T = size(IC["S"])
    end

    if ndims(extdata) != ndims(IC["S"])
        error("msfun_ica_meg_signalcorrestimate_trialwise - ERROR : External data array inconsistent.")
    end

    if epoching
        if size(extdata, 1) != K || size(extdata, 3) != T
            error("msfun_ica_meg_signalcorrestimate_trialwise - ERROR : External data not consistent with IC.")
        end
        S = size(extdata, 2)
    else
        if size(extdata, 2) != T
            error("msfun_ica_meg_signalcorrestimate_trialwise - ERROR : External data not consistent with IC.")
        end
        S = size(extdata, 1)
    end

    cfg["extname"] = get(cfg, "extname", ["ext$i" for i in 1:S])
    cfg["filter"] = get(cfg, "filter", true)
    cfg["Tsigcorr"] = get(cfg, "Tsigcorr", 0.1)
    cfg["Tpowcorr"] = get(cfg, "Tpowcorr", 0.2)

    if cfg["filter"]
        filt = get(cfg, "filt", Dict())
        filt["sfreq"] = get(filt, "sfreq", 1000)
        filt["win"] = get(filt, "win", "boxcar")
        filt["par"] = get(filt, "par", ["high", "low"])
        filt["freq"] = get(filt, "freq", [1, 25])
        filt["width"] = get(filt, "width", [0.5, 5])
        cfg["filt"] = filt
    end

    if epoching
        println("msfun_ica_meg_signalcorrestimate_trialwise - Baseline correcting and concatenating epochs...")
        X = zeros(numofic, K*T)
        Y = zeros(S, K*T)
        for k in 1:K
            avX = mean(IC["S"][k, :, :], dims=3)
            avY = mean(extdata[k, :, :], dims=3)
            X[:, (k-1)*T+1:k*T] = dropdims(IC["S"][k, :, :] .- avX, dims=3)
            Y[:, (k-1)*T+1:k*T] = dropdims(extdata[k, :, :] .- avY, dims=3)
        end
        IC["S"] = X
        extdata = Y
        T *= K
    end

    if cfg["filter"]
        println("msfun_ica_meg_signalcorrestimate_trialwise - Filtering ICs and external data...")
        win, F = msfun_filt_preparecosine(cfg["filt"], T, cfg["filt"]["sfreq"])
        IC["S"] = real(ifft(fft(IC["S"] .* win', 2) .* F, 2))
        extdata = real(ifft(fft(extdata .* win', 2) .* F, 2))
    end

    println("msfun_ica_meg_signalcorrestimate_trialwise - Performing correlation analysis...")
    sigrho = cor(extdata', IC["S']")
    powrho = cor.(extdata'.^2, IC["S'].^2)

    if !haskey(IC, "corr")
        IC["corr"] = Dict()
    end
    IC["corr"]["list"] = Int[]

    for k in 1:S
        extname = cfg["extname"][k]
        IC["corr"]["sig_IC_$(extname)"] = sigrho[k, :]
        IC["corr"]["pow_IC_$(extname)"] = powrho[k, :]
        append!(IC["corr"]["list"], findall(x -> abs(x) ≥ cfg["Tsigcorr"], sigrho[k, :]))
        append!(IC["corr"]["list"], findall(x -> abs(x) ≥ cfg["Tpowcorr"], powrho[k, :]))
    end

    IC["corr"]["list"] = sort(unique(IC["corr"]["list"]))

    println("msfun_ica_meg_signalcorrestimate_trialwise - Done.")
    return IC
end
