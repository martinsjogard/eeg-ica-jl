# ICA Tools for EEG/MEG – Julia Implementation

This repository contains a collection of Julia functions converted and adapted from MATLAB code for **Independent Component Analysis (ICA)** in EEG/MEG signal processing. These functions assist with estimating ICA components, analyzing their properties, correcting for artifacts, and visualizing results.

The prefix `msfun_` is used to denote modular signal processing functions.

## Overview

These functions support a full ICA pipeline including:

- Estimating ICA from preprocessed data
- Post hoc analyses (e.g., power spectrum, cumulants, correlation with external signals)
- Identifying artifactual components
- Visualizing component properties
- Removing artifactual components from data

---

## Function Index

### `msfun_eeg_ica_estimate.jl`
Estimates ICA from EEG data.

- **Inputs**: EEG data array (`data`), configuration dictionary (`cfg`)
- **Outputs**: Dictionary `IC` with keys like `A`, `S`, and optional diagnostics (e.g., kurtosis)

---

### `msfun_meg_ica.jl`
Estimates ICA from MEG data using a configuration struct.

- **Inputs**: MEG data (`data`), configuration (`cfg`)
- **Outputs**: ICA components in dictionary format (`IC`)

---

### `msfun_meg_ica_corranalysis.jl`  
Computes correlation between ICs and external signals.

- **Inputs**: ICA structure `IC`, external signal array `sig`, config `cfg`
- **Outputs**: Updated `IC` with `.corranalysis` field (including max correlation values)

---

### `msfun_meg_ica_corranalysis2.jl`
Computes trial-by-trial correlation between ICs and an external signal.

- **Inputs**: `IC`, external `sig`, and config
- **Outputs**: Updated `IC[:corrtrial]` with correlation structure

---

### `msfun_meg_ica_cumulantanalysis.jl`
Analyzes non-Gaussianity of ICs using 4th-order cumulants.

- **Inputs**: ICA structure
- **Outputs**: Updated `IC[:cumulant]` containing kurtosis and off-diagonal measures

---

### `msfun_meg_ica_estimate.jl`
Performs ICA decomposition on MEG data with multiple config options.

- **Inputs**: MEG signal, config
- **Outputs**: ICA result dictionary `IC`

---

### `msfun_meg_ica_ndof.jl`
Computes number of degrees of freedom for IC components.

- **Inputs**: ICA structure
- **Outputs**: Updated `IC[:ndof]` with scalar DOF per component

---

### `msfun_meg_ica_powerspectrum.jl`
Computes power spectral density for each IC.

- **Inputs**: `IC`, sampling frequency `sfreq`
- **Outputs**: Updated `IC[:freq]` and `IC[:powspctrm]`

---

### `msfun_meg_ica_resultfig.jl`
Generates result summary plots for ICA components.

- **Inputs**: `IC`, optional config
- **Outputs**: Visual output only (figures)

---

### `msfun_meg_ica_spectralanalysis.jl`
Fits linear or power-law models to IC power spectra and computes goodness-of-fit.

- **Inputs**: `IC`, fit config (`fit` types and frequency bands)
- **Outputs**: Updated `IC[:spectral]` with GOF and selection mask

---

### `msfun_meg_ica_viewer.jl`
**Placeholder**: Intended for interactive visualization of ICA results.

- **Inputs**: raw data, times, ICA struct, optional config
- **Outputs**: None (interactive GUI not implemented in Julia yet)

---

### `msfun_remartIC.jl`
Removes artifact-related ICs from EEG/MEG data and optionally visualizes before/after.

- **Inputs**: Original data (2D or 3D), `IC` struct, optional list of bad components and flag to plot
- **Outputs**: Cleaned data array (`data2`)

---

## Usage Notes

- All functions assume standard Julia packages such as `LinearAlgebra`, `Statistics`, and `Plots`.
- ICA results are stored in dictionaries using keys like `:A` (mixing matrix), `:S` (sources), and analysis results.
- Many functions can be optionally visualized using `Plots.jl`.
- For full pipelines, begin with estimation (`*_estimate.jl`), perform analyses, then use `msfun_remartIC.jl` to clean your data.