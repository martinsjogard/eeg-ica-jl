using FFTW, Statistics

function msfun_ica_meg_spectraldensity(IC::Dict, cfg::Dict)
    haskey(IC, "S") || error("msfun_ica_meg_spectraldensity - ERROR: IC structure missing elements... Try again.")
    S = IC["S"]

    ndims_S = ndims(S)
    if ndims_S == 2
        epoching = false
        (numofic, T) = size(S)
    elseif ndims_S == 3
        epoching = true
        (K, numofic, T) = size(S)
    else
        error("msfun_ica_meg_spectraldensity - ERROR: IC.S must be 2D or 3D array.")
    end

    sfreq = get(cfg, "sfreq", 1000)
    epoch_len = get(cfg, "epoch", T)
    overlap = get(cfg, "overlap", 2)

    println("msfun_ica_meg_spectraldensity - Computing Fourier power spectrum of ICs...")
    println("msfun_ica_meg_spectraldensity -       epoch length $(epoch_len/sfreq) sec...")
    println("msfun_ica_meg_spectraldensity -       epochs overlap $overlap...")

    if !epoching
        powspctrm = zeros(numofic, epoch_len)
        ti, tf = 1, epoch_len
        nav = 0
        while tf <= T
            epoch = S[:, ti:tf]
            powspctrm .+= abs.(fft(epoch, 2)).^2
            nav += 1
            ti += round(Int, epoch_len / overlap)
            tf = ti + epoch_len - 1
        end
        powspctrm ./= nav
    else
        powspctrm = zeros(K, numofic, epoch_len)
        ti, tf = 1, epoch_len
        nav = 0
        while tf <= T
            epoch = S[:, :, ti:tf]
            powspctrm .+= abs.(fft(epoch, dims=3)).^2
            nav += 1
            ti += round(Int, epoch_len / overlap)
            tf = ti + epoch_len - 1
        end
        powspctrm = dropdims(mean(powspctrm ./ nav, dims=1), dims=1)
    end

    freq = sfreq * collect(0:floor(Int, epoch_len/2)-1) ./ epoch_len
    IC["powspctrm"] = powspctrm[:, 1:floor(Int, epoch_len/2)] ./ epoch_len
    IC["freq"] = freq

    println("msfun_ica_meg_spectraldensity -       frequency domain [$(freq[1]) $(freq[end])] Hz...")
    println("msfun_ica_meg_spectraldensity - Done.")

    return IC
end
