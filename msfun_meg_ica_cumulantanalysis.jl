using Statistics

function msfun_meg_ica_cumulantanalysis(IC::Dict, cfg::Dict)
    if !haskey(IC, "S") || !(IC["S"] isa AbstractArray) || !(ndims(IC["S"]) in (2, 3))
        error("msfun_meg_ica_cumulantanalysis - ERROR : IC structure missing elements or inconsistent.")
    end

    epoching = ndims(IC["S"]) == 3
    if epoching
        K, numofic, T = size(IC["S"])
    else
        numofic, T = size(IC["S"])
    end

    Tskew = get(cfg, "Tskew", NaN)
    Tkurt = get(cfg, "Tkurt", 15.0)

    if epoching
        println("msfun_meg_ica_cumulantanalysis - Baseline correcting and concatenating epochs...")
        X = zeros(numofic, K*T)
        for k in 1:K
            av = mean(IC["S"][k, :, :], dims=3)
            X[:, (k-1)*T+1:k*T] .= dropdims(IC["S"][k, :, :] .- av, dims=3)
        end
        IC["S"] = X
    end

    println("msfun_meg_ica_cumulantanalysis - Computing IC skewness and kurtosis...")
    skew_vals = [Statistics.skewness(IC["S"][i, :]) for i in 1:numofic]
    kurt_vals = [Statistics.kurtosis(IC["S"][i, :]) for i in 1:numofic]

    println("msfun_meg_ica_cumulantanalysis - Classifying ICs using their kurtosis...")
    sorted_indices = sortperm(kurt_vals, rev=true)
    IC["S"] = IC["S"][sorted_indices, :]
    if haskey(IC, "A")
        IC["A"] = IC["A"][:, sorted_indices]
    end
    if haskey(IC, "W")
        IC["W"] = IC["W"][sorted_indices, :]
    end

    cumulant = Dict(
        "skew" => skew_vals[sorted_indices],
        "kurt" => kurt_vals[sorted_indices],
        "Tskew" => Tskew,
        "Tkurt" => Tkurt
    )

    list_skew = isnan(Tskew) ? Int[] : findall(abs.(cumulant["skew"]) .> Tskew)
    list_kurt = isnan(Tkurt) ? Int[] : findall(cumulant["kurt"] .> Tkurt)
    cumulant["list"] = sort(unique(vcat(list_skew, list_kurt)))

    IC["cumulant"] = cumulant
    println("msfun_meg_ica_cumulantanalysis - Done.")
    return IC
end