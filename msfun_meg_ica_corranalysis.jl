using FFTW, Statistics, LinearAlgebra

include("msfun_prepare_cosine_filter.jl")

function msfun_meg_ica_corranalysis(IC::Dict, extdata, cfg::Dict)
    if extdata === nothing || isempty(extdata)
        println("msfun_meg_ica_corranalysis - WARNING : No external data supplied... Skipping the correlation analysis.")
        return IC
    end

    S_data = IC["S"]
    if ndims(S_data) == 3
        epoching = true
        K, numofic, T = size(S_data)
    else
        epoching = false
        numofic, T = size(S_data)
    end

    extdata = Array(extdata)
    if epoching && size(extdata, 1) != K
        error("Epoch count mismatch between IC and extdata")
    end

    S = epoching ? size(extdata, 2) : size(extdata, 1)
    cfg["extname"] = get(cfg, "extname", ["ext\$k" for k in 1:S])
    cfg["filter"] = get(cfg, "filter", true)
    cfg["Tcorr"] = get(cfg, "Tcorr", 0.15)

    if cfg["filter"]
        cfg["filt"] = get(cfg, "filt", Dict(
            "sfreq" => 1000,
            "win" => "boxcar",
            "par" => ["high", "low"],
            "freq" => [1.0, 25.0],
            "width" => [0.5, 5.0]
        ))
    end

    if epoching
        println("msfun_meg_ica_corranalysis - Baseline correcting and concatenating epochs...")
        X = zeros(numofic, K * T)
        Y = zeros(S, K * T)
        for k in 1:K
            avX = dropdims(mean(S_data[k, :, :], dims=2), dims=2)
            avY = dropdims(mean(extdata[k, :, :], dims=2), dims=2)
            X[:, (k-1)*T+1:k*T] = S_data[k, :, :] .- avX
            Y[:, (k-1)*T+1:k*T] = extdata[k, :, :] .- avY
        end
        IC["S"] = X
        extdata = Y
    end

    if cfg["filter"]
        println("msfun_meg_ica_corranalysis - Filtering ICs and external data...")
        win, F = msfun_prepare_cosine_filter(cfg["filt"], size(IC["S"], 2), cfg["filt"]["sfreq"])
        icasig = real(ifft(fft(IC["S"] .* win, 2) .* F, 2))
        extdata = real(ifft(fft(extdata .* win, 2) .* F, 2))
    else
        icasig = IC["S"]
    end

    println("msfun_meg_ica_corranalysis - Performing correlation analysis...")
    rho = cor(transpose(extdata), transpose(icasig))

    IC["corr"] = Dict{String, Any}()
    IC["corr"]["list"] = Int[]

    for k in 1:S
        extname = cfg["extname"][k]
        IC["corr"][extname] = rho[k, :]
        IC["corr"]["list"] = union(IC["corr"]["list"], findall(x -> abs(x) >= cfg["Tcorr"], rho[k, :]))
    end

    IC["corr"]["Tcorr"] = cfg["Tcorr"]
    println("msfun_meg_ica_corranalysis - Done.")
    return IC
end