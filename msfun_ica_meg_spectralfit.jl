using LinearAlgebra
using Plots
using Statistics
using Printf

function msfun_ica_meg_spectralfit(IC::Dict, cfg::Dict)
    if !haskey(IC, "freq") || !haskey(IC, "powspctrm")
        error("IC structure missing required fields 'freq' and 'powspctrm'")
    end

    powspctrm = IC["powspctrm"]
    numofic, F = size(powspctrm)
    freq = IC["freq"]
    IC["spectral"] = Dict("fit" => copy(cfg["fit"]))
    L = Int(length(cfg["fit"]) ÷ 2)
    gof = zeros(numofic, L)

    for k in 1:L
        kind = lowercase(cfg["fit"][2k - 1])
        fmin, fmax = cfg["fit"][2k]
        f1 = findfirst(x -> x >= fmin, freq)
        f2 = findlast(x -> x <= fmax, freq)
        cfg["fit"][2k] = [f1, f2]
        IC["spectral"]["fit"][2k] = freq[f1:f2]
    end

    if !haskey(cfg, "Tgof")
        cfg["Tgof"] = fill(0.03, L)
    end

    visual = get(cfg, "visual", true)
    if visual
        plt = plot(layout = (2, L), size=(400*L, 600))
    end

    for n in 1:numofic
        for k in 1:L
            fidx = cfg["fit"][2k]
            x = powspctrm[n, fidx[1]:fidx[2]]
            f = freq[fidx[1]:fidx[2]]
            xhat = similar(x)

            if lowercase(cfg["fit"][2k - 1]) == "linear"
                X = hcat(ones(length(f)), f)
                β = X \ x
                xhat = X * β
            elseif lowercase(cfg["fit"][2k - 1]) == "powlaw"
                logf = log.(f)
                logx = log.(x)
                coef = polyfit(logf, logx, 1)
                a_init = exp(coef[2])
                b_init = coef[1]
                xhat = a_init .* f .^ b_init
            else
                error("Unsupported fit type: $(cfg["fit"][2k - 1])")
            end

            gof[n, k] = sum((x .- xhat).^2) / sum(x.^2)

            if visual
                plot!(plt[1, k], f, x, label="data")
                plot!(plt[1, k], f, xhat, label="fit", color=:red)
                title!(plt[1, k], "IC $n - $(cfg["fit"][2k - 1]) fit")
                err = (x .- xhat).^2 ./ mean(x.^2)
                scatter!(plt[2, k], f, err, label="", color=:red)
                hline!(plt[2, k], [cfg["Tgof"][k]], linestyle=:dash, label="")
                title!(plt[2, k], @sprintf("error (mean = %.3f)", gof[n, k]))
            end
        end
    end

    IC["spectral"]["gof"] = gof
    IC["spectral"]["Tgof"] = cfg["Tgof"]
    IC["spectral"]["list"] = [i for i in 1:numofic if all(gof[i, j] < cfg["Tgof"][j] for j in 1:L)]

    if visual
        display(plt)
    end

    return IC
end
