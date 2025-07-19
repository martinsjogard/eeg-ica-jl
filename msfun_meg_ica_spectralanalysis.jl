
using LinearAlgebra, Statistics, Printf
using Plots

function msfun_meg_ica_spectralanalysis(IC::Dict, cfg::Dict)
    freq = vec(IC["freq"])
    powspctrm = IC["powspctrm"]
    numofic, F = size(powspctrm)

    if size(powspctrm, 2) != length(freq)
        error("Frequency and power spectrum dimensions do not match")
    end

    cfg = copy(cfg)
    cfg["visual"] = get(cfg, "visual", true)
    fitcfg = cfg["fit"]
    L = length(fitcfg) ÷ 2
    cfg["Tgof"] = get(cfg, "Tgof", fill(0.03, L))

    IC["spectral"] = Dict()
    IC["spectral"]["fit"] = deepcopy(fitcfg)

    for k in 1:L
        fband = fitcfg[2k]
        idx_start = findfirst(>=(fband[1]), freq)
        idx_end = findlast(<=(fband[2]), freq)
        IC["spectral"]["fit"][2k] = [freq[idx_start], freq[idx_end]]
    end

    IC["spectral"]["gof"] = zeros(numofic, L)

    if cfg["visual"]
        display(Plots.plot(layout=(2, L)))
    end

    for n in 1:numofic
        for k in 1:L
            fit_type = fitcfg[2k - 1]
            fband = fitcfg[2k]
            idx_start = findfirst(>=(fband[1]), freq)
            idx_end = findlast(<=(fband[2]), freq)
            x = powspctrm[n, idx_start:idx_end]
            f = freq[idx_start:idx_end]

            if fit_type == "linear"
                p = polyfit(f, x, 1)
                xhat = p[1] .* f .+ p[2]
            elseif fit_type == "powlaw"
                p = polyfit(log.(f), log.(x), 1)
                xhat = exp(p[2]) .* f .^ p[1]
            else
                error("Unsupported fit type")
            end

            IC["spectral"]["gof"][n, k] = sum((x .- xhat).^2) / sum(x.^2)

            if cfg["visual"]
                plot!(f, x, label="IC $n", subplot=k)
                plot!(f, xhat, color=:red, subplot=k)
                plot!(f, (x .- xhat).^2 ./ mean(x.^2), seriestype=:scatter, subplot=L+k)
                hline!([cfg["Tgof"][k]], linetype=:dash, color=:black, subplot=L+k)
            end
        end
    end

    mask = trues(numofic)
    for k in 1:L
        mask .&= IC["spectral"]["gof"][:, k] .< cfg["Tgof"][k]
    end
    IC["spectral"]["list"] = findall(mask)
    IC["spectral"]["Tgof"] = cfg["Tgof"]

    return IC
end
