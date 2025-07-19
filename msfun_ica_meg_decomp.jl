using LinearAlgebra, Statistics, MultivariateStats

function msfun_ica_meg_decomp(sig::Array{Float64,2}, cfg::Dict=Dict())
    N, T = size(sig)

    normalize = get(cfg, "normalize", nothing)
    if normalize !== nothing
        normalize = collect(normalize)
        if length(normalize) != N || any(x -> x <= 0, normalize)
            error("msfun_ica_meg_decomp - ERROR: Normalization factors inconsistent... Try again.")
        end
    else
        println("msfun_ica_meg_decomp - using temporal standard deviation of data...")
        normalize = std(sig, dims=2)[:]
    end

    println("msfun_ica_meg_decomp - Normalizing data units...")
    sig_norm = sig ./ normalize

    println("msfun_ica_meg_decomp - FASTICA in action...")
    ica_model = fastica(transpose(sig_norm), maxiter=200, fun="tanh")
    S = transpose(ica_model.S)
    A = ica_model.A
    W = ica_model.W

    println("msfun_ica_meg_decomp - Restoring original data units...")
    A = A .* normalize
    W = W ./ transpose(normalize)

    IC = Dict(
        "S" => S,
        "A" => A,
        "W" => W
    )

    println("msfun_ica_meg_decomp - Done.")
    return IC
end
