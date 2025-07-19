
using LinearAlgebra, Statistics, Printf, Plots

function msfun_meg_ica_ndof(data::Array, cfg::Dict)
    ndims(data) in [2, 3] || error("msfun_meg_ica_ndof - ERROR: data must be a numeric array... Try again.")

    epoching = ndims(data) == 3
    if epoching
        K, N, T = size(data)
    else
        N, T = size(data)
    end

    normalize = get(cfg, "normalize", nothing)
    if normalize !== nothing
        normalize = collect(normalize)
        length(normalize) == N || error("msfun_meg_ica_ndof - ERROR: Normalization factors inconsistent... Try again.")
        any(x -> x <= 0, normalize) && error("msfun_meg_ica_ndof - ERROR: Normalization values must be positive.")
    end

    method = get(cfg, "method", "rel")
    method in ["abs", "maxrel", "rel"] || error("msfun_meg_ica_ndof - ERROR: Method not recognized... Try again.")

    param = get(cfg, "param", 1e3)

    if epoching
        println("msfun_meg_ica_ndof - Baseline correcting and concatenating epochs...")
        X = zeros(N, K*T)
        for k in 1:K
            av = mean(data[k, :, :], dims=3)[:]
            X[:, (k-1)*T+1:k*T] = dropdims(data[k, :, :], dims=1) .- av
        end
        data = X
    end

    println("msfun_meg_ica_ndof - Normalizing data...")
    if normalize === nothing
        normalize = std(data, dims=2)[:]
        println("msfun_meg_ica_ndof -   using data standard deviation...")
    end
    data = data ./ normalize

    println("msfun_meg_ica_ndof - Computing normalized data covariance and its eigenvalues...")
    D = eigvals(cov(data'))
    D = sort(D)

    if method == "abs"
        cutoff = param
        n = findfirst(x -> x >= cutoff, D)
    elseif method == "maxrel"
        cutoff = maximum(D) / param
        n = findfirst(x -> x >= cutoff, D)
    else
        R = D[2:end] ./ D[1:end-1]
        n = findlast(x -> x >= param, R)
    end

    ndof = N - n
    println("msfun_meg_ica_ndof - Estimated $ndof largest eigendirections...")

    if method in ["abs", "maxrel"]
        plot(D, label="Eigenvalues", color=:red)
        plot!((n+1):N, D[n+1:end], label="Below cutoff")
        hline!([cutoff], color=:green, label="Cutoff")
    else
        plot(layout=(1, 2))
        plot!(1, D, label="Eigenvalues", xlabel="n", ylabel="eigenvalue(n)", color=:red)
        plot!(1, (n+1):N, D[n+1:end], label="Below cutoff")
        plot!(2, R, label="Ratios", xlabel="n", ylabel="eigenvalue(n+1)/eigenvalue(n)", color=:red)
        plot!(2, (n+1):(length(R)), R[n+1:end], label="Below", color=:blue)
        hline!([param], subplot=2, color=:green)
    end

    go = readline(stdin, keep=true, prompt="Keep $ndof largest eigendirections? [y/n] ")
    if lowercase(strip(go)) in ["n", "no"]
        print("Enter number of eigendirections to keep: ")
        ndof = parse(Int, readline())
    end

    return ndof
end
