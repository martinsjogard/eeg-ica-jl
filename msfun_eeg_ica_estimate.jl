using LinearAlgebra
using Statistics
using FastICA  # Assumes ICA.jl or equivalent package is installed

function msfun_eeg_ica_estimate(sig::Matrix{Float64}, cfg::Dict=Dict())
    N, T = size(sig)
    if ndims(sig) != 2
        error("sig must be a 2D array [N x T]")
    end

    normalize = get(cfg, "normalize", nothing)
    fastica_params = get(cfg, "fastica", Dict("fun" => "tanh", "n_components" => min(N, 30)))

    println("msfun_eeg_ica_estimate - Normalizing data units...")
    if normalize === nothing
        println("msfun_eeg_ica_estimate -       using temporal standard deviation of data...")
        normalize = std(sig, dims=2)
    else
        normalize = vec(normalize)
        if length(normalize) != N || any(normalize .<= 0)
            error("normalize must be a positive vector of length N")
        end
    end

    sig_norm = sig ./ normalize

    println("msfun_eeg_ica_estimate - FASTICA in action...")
    num_components = fastica_params["n_components"]
    fun_type = fastica_params["fun"]

    result = fastica(sig_norm; n_components=num_components, fun=fun_type)
    S = result.S
    A = result.A
    W = result.W

    println("msfun_eeg_ica_estimate - Restoring original data units...")
    A_scaled = A .* normalize
    W_scaled = W ./ normalize'

    println("msfun_eeg_ica_estimate - Done.")

    return Dict("A" => A_scaled, "W" => W_scaled, "S" => S)
end