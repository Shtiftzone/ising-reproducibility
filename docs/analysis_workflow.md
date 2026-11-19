# Analysis and figure reproduction

This document describes the statistical-analysis and plotting workflow for the two-dimensional Ising model. It covers the transformation of raw simulation measurements and multihistogram results into processed CSV files, finite-size-scaling (FSS) estimates, and publication figures.

The simulation and multihistogram stages are documented separately in [simulation_workflow.md](simulation_workflow.md) and [multihistogram_workflow.md](multihistogram_workflow.md).

## Scope and directory conventions

The analysis accompanying the paper uses two lattice/representation combinations:

- `square/cell`
- `triangular/vertex`

Scripts may support other combinations, but the commands below reproduce the analysis for these two datasets.

Run all commands from the repository root. Python dependencies are listed in `requirements.txt` and can be installed with:

```bash
python3 -m pip install -r requirements.txt
```

The principal directories are:

```text
results/square_cell_simulations/                 Raw square/cell measurements (external)
results/triangular_vertex_simulations/          Raw triangular/vertex measurements (external)
results/multihistogram/<lattice>/<representation>/
                                                Multihistogram outputs (external)
data/analysis_csv/<lattice>/<representation>/  Processed analysis CSVs
data/table_csv/<lattice>/<representation>/     FSS tables
results/figures/<lattice>/<representation>/   Generated figures
```

Large raw simulation time series and multihistogram outputs are distributed separately. Selected processed CSVs and figures are included in the repository, so the plotting scripts can be run without repeating the production simulations.

## Workflow overview

```text
Raw simulation measurements                  Multihistogram bootstrap results
           |                                             |
           +--> Temperature-dependent                Pseudocritical peaks
           |    observables and block jackknife           |
           |                                             +--> Leading FSS scans
           +--> Distribution moments                     |    and stability plots
           |    and effective resolution                  |
           |                                             +--> Square-lattice
           |                                                  subleading FSS
           v                                             v
    Processed analysis CSVs --------------------> Figures and tables
```

All processed data are labeled by lattice type and representation. When changing datasets, specify the corresponding input directory explicitly; the labels alone do not verify the physical origin of an input file.

## 1. Temperature-dependent observables

**Analysis script:** `scripts/analysis/compute_observables_vs_temperature.py`

The script reads the recorded measurements for a selected lattice size from:

```text
<results_dir>/T_<temperature>/size_<L>/mefile.txt
<results_dir>/T_<temperature>/size_<L>/epfile.txt
<results_dir>/T_<temperature>/size_<L>/enfile.txt
```

It calculates the temperature dependence of:

- Mean absolute magnetization, `<|M|>`.
- Mean energy, `<H>`.
- Mean absolute symmetrized Euler characteristic, `<|EC_sym|>`.
- Mean averaged Euler characteristic, `<EC_avg>`.

Here `EC_sym = (EC_minus - EC_plus)/2` and `EC_avg = (EC_minus + EC_plus)/2`. The sign convention for `EC_sym` does not affect its absolute value.

### Statistical uncertainties

The script estimates the integrated autocorrelation time separately for the four measured observable series at each temperature. It selects a common block length from the largest estimate and computes a **nonoverlapping, delete-one-block jackknife** standard error for each mean.

The example below uses `--tau-multiplier 10 --min-blocks 10`, as in the analysis workflow. If there are insufficient measurements to attain the target block length while preserving the minimum number of blocks, the block length is capped and a warning is issued. Check `autocorrelation_diagnostics.csv` for the estimated autocorrelation times, selected block sizes, discarded trailing samples, and comparisons with unblocked standard errors.

These are **block-jackknife** errors, not the paired bootstrap errors used for subleading FSS.

### Commands

Square lattice, `L=3072`:

```bash
python3 scripts/analysis/compute_observables_vs_temperature.py \
    --results-dir results/square_cell_simulations \
    --lattice-type square --representation cell \
    --size 3072 --output-dir data/analysis_csv \
    --tau-multiplier 10 --min-blocks 10
```

Triangular lattice, `L=3072`:

```bash
python3 scripts/analysis/compute_observables_vs_temperature.py \
    --results-dir results/triangular_vertex_simulations \
    --lattice-type triangular --representation vertex \
    --size 3072 --output-dir data/analysis_csv \
    --tau-multiplier 10 --min-blocks 10
```

The output is stored in:

```text
data/analysis_csv/<lattice>/<representation>/observables_vs_temperature/size_3072/
    mag_vs_T.csv
    energy_vs_T.csv
    euler_sym_vs_T.csv
    ec_avg_vs_T.csv
    autocorrelation_diagnostics.csv
```

Each observable CSV contains the columns `L`, `T`, `observable`, `mean`, `jackknife_se`, `n_samples`, `n_used`, `block_size`, `n_blocks`, and `tau_int`, in addition to metadata. Use `--size` to analyze another lattice size.

The corresponding plotting script is `scripts/plotting/plot_observables_vs_temperature.py` (see Section 6).

## 2. Pseudocritical temperatures from multihistogram results

**Analysis script:** `scripts/analysis/compute_pseudocritical_temperatures.py`

The script reads:

```text
results/multihistogram/<lattice>/<representation>/results_<L>.dat
```

For each lattice size and bootstrap ID, it finds the maximum of the susceptibility-like quantity for three observables:

- `ec`: symmetrized Euler-characteristic observable.
- `m`: magnetization observable.
- `d`: related observable used in the paper.

The extraction reads the multihistogram temperature from column 2 and divides the corresponding stored susceptibility-related column by `T` before locating its maximum. The zero-based input column indices for `ec`, `m`, and `d` are 4, 6, and 8, respectively.

For each observable and bootstrap realization, the peak location defines `Tc`; the peak height is recorded as `chi_c`. The bootstrap IDs are preserved so the observables can be compared using paired samples later.

### Commands

```bash
for combo in "square cell" "triangular vertex"; do
    read -r LAT REP <<< "$combo"
    python3 scripts/analysis/compute_pseudocritical_temperatures.py \
        --input-dir "results/multihistogram/$LAT/$REP" \
        --lattice-type "$LAT" --representation "$REP" \
        --output-dir data/analysis_csv
done
```

This produces six CSVs per combination:

```text
data/analysis_csv/<lattice>/<representation>/
    pseudocritical_ec.csv
    pseudocritical_m.csv
    pseudocritical_d.csv
    pseudocritical_bootstrap_ec.csv
    pseudocritical_bootstrap_m.csv
    pseudocritical_bootstrap_d.csv
```

The summary CSVs include `L`, the bootstrap mean `Tc`, its sample standard deviation `Tc_std`, the mean peak height `chi_c`, its standard deviation `chi_c_std`, and the number of bootstrap realizations `N_boot`.

**Important:** `Tc_std` is the sample standard deviation of bootstrap peak locations, not the standard error of their bootstrap mean. The FSS scripts use `Tc_std` as the uncertainty for each lattice size.

## 3. Leading finite-size scaling

The leading FSS model is

$$
T_c(L) = T_c(\infty) + A L^{-1/\nu}.
$$

Each observable (`ec`, `m`, `d`) is fitted separately, with `Tc_inf`, `A`, and `nu` free. The fits use the mean `Tc(L)` and its `Tc_std` as inputs, weighted using the supplied standard deviations (`absolute_sigma=True` in SciPy).

Fit quality is assessed using the chi-square statistic, `chi2/dof`, and a chi-square goodness-of-fit p-value. **AIC and AICc are not used in the current workflow.**

### Scan over minimum lattice size

**Script:** `scripts/analysis/compute_fss_lmin.py`

The script repeats the leading FSS fit for successive lower cutoffs `L_min`, while retaining the available upper end of the data.

```bash
for combo in "square cell" "triangular vertex"; do
    read -r LAT REP <<< "$combo"
    for OBS in ec m d; do
        python3 scripts/analysis/compute_fss_lmin.py \
            --input "data/analysis_csv/$LAT/$REP/pseudocritical_${OBS}.csv" \
            --output "data/table_csv/$LAT/$REP/lmin_${OBS}.csv"
    done
done
```

Results are written to `data/table_csv/<lattice>/<representation>/lmin_{ec,m,d}.csv`. The tables record fitted parameters and uncertainties together with goodness-of-fit diagnostics.

### Stability scans at two fixed upper cutoffs

**Script:** `scripts/analysis/compute_fss_fixed_lmax.py`

This script repeats the fits while varying `L_min` for the largest and second-largest lattice sizes available in the common dataset. With the full 12-size dataset these upper cutoffs are `L_max=3072` and `L_max=2048`. The default minimum number of fitted sizes is six.

```bash
for combo in "square cell" "triangular vertex"; do
    read -r LAT REP <<< "$combo"
    python3 scripts/analysis/compute_fss_fixed_lmax.py \
        --lattice-type "$LAT" --representation "$REP" \
        --input-dir data/analysis_csv \
        --output-dir data/table_csv \
        --min-points 6
done
```

The files are saved to:

```text
data/table_csv/<lattice>/<representation>/fixed_lmax/
    lmin_ec_lmax_max.csv
    lmin_m_lmax_max.csv
    lmin_d_lmax_max.csv
    lmin_ec_lmax_second.csv
    lmin_m_lmax_second.csv
    lmin_d_lmax_second.csv
```

The plotting script `plot_fss_lmin_stability.py` reads these CSVs to produce four panels per lattice/representation combination: `Tc_inf` and `nu` versus `L_min`, each for both upper cutoffs.

## 4. Subleading FSS corrections on the square lattice

**Script:** `scripts/analysis/paired_fss_corrections.py`

For `square/cell` only, the analysis also uses the model with fixed `nu=1`:

$$
T_c(L) = T_c(\infty) + \frac{a}{L} + \frac{b}{L^2}.
$$

The central estimates of `Tc_inf`, `a`, and `b` are obtained by weighted linear least squares using the **mean** pseudocritical temperature at each `L`, with weights proportional to `1/Tc_std^2`.

The uncertainty analysis draws from the *existing* bootstrap peak-location samples. At each lattice size, `ec`, `m`, and `d` are sampled jointly using their shared bootstrap IDs; different lattice sizes are sampled independently. The three-parameter fit is repeated for each paired draw to estimate parameter uncertainties and the distributions of differences such as `b_ec - b_m` and `b_ec - b_d`. This is a secondary paired resampling of extracted peaks, **not** a new moving-block bootstrap of the raw Monte Carlo time series.

```bash
python3 scripts/analysis/paired_fss_corrections.py \
    --lattice-type square --representation cell \
    --input-base data/analysis_csv \
    --output-base data/analysis_csv/fss_corrections \
    --lmin 64 --lmax 3072 \
    --n-resamples 10000 --seed 123456 --ci 95
```

Outputs:

```text
data/analysis_csv/fss_corrections/square/cell/
    fss_subleading_fit_table.csv
    paired_fss_bootstrap.csv
    paired_fss_summary.csv
    delta_b_ec_m.png
    delta_b_ec_d.png
    delta_abs_b_ec_m.png
    delta_abs_b_ec_d.png
```

`fss_subleading_fit_table.csv` reports the central parameters, bootstrap standard deviations and confidence intervals, weighted-least-squares uncertainties, chi-square diagnostics, and p-values. `paired_fss_summary.csv` contains the bootstrap distributions of `b` and its paired differences, including percentile confidence intervals. The four PNGs visualize those comparison distributions.

The subleading correction fit is not part of the triangular-lattice workflow used in the paper.

## 5. Distribution shape and effective resolution

**Analysis script:** `scripts/analysis/compute_distribution_resolution.py`

For a selected **measured** temperature, the script reads aligned raw measurements of `M`, `EC_plus`, and `EC_minus` from `mefile.txt`, `epfile.txt`, and `enfile.txt`, respectively. It constructs

$$
EC_{\mathrm{sym}} = \frac{EC_{\mathrm{minus}} - EC_{\mathrm{plus}}}{2},
$$

and

$$
D = \begin{cases}
M + 2 EC_{\mathrm{sym}}, & \text{square lattice}, \\
M + EC_{\mathrm{sym}}, & \text{triangular lattice}.
\end{cases}
$$

The script writes two datasets:

- `distribution_shape.csv`: sample means, standard deviations, skewness, and excess kurtosis for `M`, `EC_sym`, and `D`.
- `entropy_effective_resolution.csv`: observed unique-value counts, Shannon entropies, entropy-based effective numbers of values, and ratios relative to magnetization.

For empirical probabilities `p_i`, Shannon entropy and effective resolution are

$$
H = -\sum_i p_i\log p_i, \qquad K_{\mathrm{eff}} = e^H.
$$

The unique-value count is denoted by `K`. The plotted ratios are `K_EC/K_M`, `K_D/K_M`, `Keff_EC/Keff_M`, and `Keff_D/Keff_M`.

**Kolmogorov-Smirnov tests and KS plots have been removed from the current workflow.** Distribution moments are retained in `distribution_shape.csv` for descriptive analysis.

### Commands

The processed data currently tracked in the repository use `T=2.26900` for `square/cell` and `T=3.64100` for `triangular/vertex`.

```bash
python3 scripts/analysis/compute_distribution_resolution.py \
    --results-dir results/square_cell_simulations \
    --lattice-type square --representation cell \
    --temperature 2.26900 --output-dir data/analysis_csv

python3 scripts/analysis/compute_distribution_resolution.py \
    --results-dir results/triangular_vertex_simulations \
    --lattice-type triangular --representation vertex \
    --temperature 3.64100 --output-dir data/analysis_csv
```

Output files are placed in:

```text
data/analysis_csv/<lattice>/<representation>/distribution_resolution/
    distribution_shape.csv
    entropy_effective_resolution.csv
```

The plotting script reads only `entropy_effective_resolution.csv` and does not generate KS-distance figures.

## 6. Figure generation from processed CSVs

**All commands below use existing processed CSV files.** They do not rerun the simulations, multihistogram reweighting, or statistical-analysis scripts.

```bash
for combo in "square cell" "triangular vertex"; do
    read -r LAT REP <<< "$combo"

    python3 scripts/plotting/plot_distribution_resolution.py \
        --lattice-type "$LAT" --representation "$REP"

    python3 scripts/plotting/plot_fss_lmin_stability.py \
        --lattice-type "$LAT" --representation "$REP"

    python3 scripts/plotting/plot_fss_tc.py \
        --lattice-type "$LAT" --representation "$REP"

    python3 scripts/plotting/plot_observables_vs_temperature.py \
        --lattice-type "$LAT" --representation "$REP" --size 3072
done
```

The outputs have the following structure:

```text
results/figures/<lattice>/<representation>/
    <lattice>_<representation>_resolution_ratios_vs_L.png
    <lattice>_<representation>_fss_Tc_ec_sym_vs_m_vs_d.png
    Tc_vs_Lmin_<lattice>_<representation>_Lmax_max.png
    Tc_vs_Lmin_<lattice>_<representation>_Lmax_second.png
    nu_vs_Lmin_<lattice>_<representation>_Lmax_max.png
    nu_vs_Lmin_<lattice>_<representation>_Lmax_second.png
    observables_vs_temperature/size_3072/
        mag_vs_T_abs.png
        energy_vs_T.png
        euler_sym_vs_T_abs.png
        ec_avg_vs_T.png
```

`plot_fss_tc.py` performs its own leading FSS fits from the pseudocritical summary CSVs. Its per-observable `L_min` cutoffs can be overridden with `--lmin-ec`, `--lmin-m`, and `--lmin-d`. For the current paper combinations, the script defaults to `64/64/64` for `square/cell` and `96/64/64` for `triangular/vertex` (ordered `ec/m/d`). The default horizontal coordinate is `1/L`.

The theoretical values used as reference lines on the stability plots are

$$
\nu = 1, \quad T_c^{\mathrm{square}} = \frac{2}{\ln(1+\sqrt{2})},
\quad T_c^{\mathrm{triangular}} = \frac{4}{\ln 3}.
$$

Figure files are now written directly to their lattice/representation subdirectories. A separate `organize_figures.sh` step is not required.

## 7. Statistical diagnostics

The following scripts provide additional checks of temporal autocorrelation, proposed block lengths, correlations between simulation streams, and the overlap of adjacent energy histograms. These diagnostics are separate from the production multihistogram and FSS calculations. They do not change the published estimates or rerun the Monte Carlo simulations.

Two different input layouts are used:

```text
# Prepared multihistogram input (ten whitespace-separated columns):
data/eul2d/<lattice>/<representation>/size_<L>/conf-<T>-<run>.dat

# Raw simulation output:
results/<simulation_type>/T_<T>/size_<L>/mefile.txt
```

The `conf-<T>-<run>.dat` files are generated by `scripts/multihistogram/prepare_eul2d_input.sh` from the original measurements. Their columns begin with magnetization, energy, positive Euler characteristic, and negative Euler characteristic, followed by six cell-count columns. The autocorrelation and cross-stream scripts below use these prepared files; the energy-overlap script reads `mefile.txt` directly. To analyze a different combination, replace `square/cell` with the appropriate lattice and representation, and update the `--lattice-type` or `--lattice` argument.

### 7.1. Integrated autocorrelation times

**Script:** `scripts/analysis/diagnose_autocorr.py`

This script estimates integrated autocorrelation times for seven recorded series:

- Energy per spin.
- Absolute magnetization per spin.
- Absolute symmetrized Euler characteristic per spin.
- Absolute domain-imbalance observable `D` per spin.
- The squares of the normalized magnetization, symmetrized Euler characteristic, and `D`.

The signed observable used to construct `D` is `D = M + 2 EC_sym` on the square lattice and `D = M + EC_sym` on the triangular lattice, with `EC_sym = (EC_plus - EC_minus)/2` in these diagnostic scripts. The sign convention for `EC_sym` is relevant when combining it with `M`.

The script compiles a temporary C++ driver using the repository's original `src/multihistogram/Autocorr.cpp` and `Autocorr.h`. It analyzes each temperature, size, and run separately; recorded time series are not concatenated.

For example:

```bash
python3 scripts/analysis/diagnose_autocorr.py \
    --input-dir data/eul2d/square/cell \
    --autocorr-source src/multihistogram/Autocorr.cpp \
    --lattice-type square \
    --output-dir data/analysis_csv/diagnostics/square/cell/autocorr
```

The principal outputs are:

```text
autocorr_all.csv
autocorr_summary_by_L.csv
autocorr_invalid.csv
autocorr_summary.txt
autocorr_vs_T_<observable>.png
```

In the estimator convention used here, the independent-sample reference is `tau_int = 0.5`, with an approximate effective sample count `N_eff = N/(2*tau_int)`. Autocorrelation times are measured in **intervals between saved measurements**, not individual Metropolis sweeps or Wolff updates. These estimates are diagnostic: short series and the autocorrelation window can miss slow modes, and these results do not by themselves validate uncertainties of multihistogram peak estimates.

### 7.2. Block-length feasibility for block bootstrap

**Script:** `scripts/analysis/analyze_block_bootstrap.py`

This extends the autocorrelation checks above by proposing a **common block length for each time series**, based on the largest valid `tau_int` among the seven observables. With block factor `f`, the proposed nonoverlapping block length is

```text
block_length = ceil(f * max(0.5, tau_max))
n_blocks     = floor(N / block_length)
```

The default factor is `10`, and the default minimum number of blocks is `20`. The following example explicitly requests a minimum of `10` blocks:

```bash
python3 scripts/analysis/analyze_block_bootstrap.py \
    --input-dir data/eul2d/square/cell \
    --autocorr-source src/multihistogram/Autocorr.cpp \
    --lattice-type square \
    --output-dir data/analysis_csv/diagnostics/square/cell/block_feasibility \
    --block-factor 10 --min-blocks 10
```

The outputs include:

```text
autocorr_all.csv
autocorr_invalid.csv
autocorr_summary_by_L.csv
block_sizes_by_series.csv
block_problem_cases.csv
block_summary_by_L.csv
block_summary.txt
autocorr_vs_T_<observable>.png
tau_max_vs_L.png
block_length_vs_L.png
n_blocks_vs_L.png
n_blocks_histogram.png
```

The status fields distinguish feasible series from cases with too few blocks or incomplete autocorrelation estimates. The `tau_int >= 0.5` floor is used **only to propose block lengths**; the raw estimates are preserved. **This script does not perform bootstrap resampling or multihistogram reweighting.** Proposed block lengths should be checked for statistical stability before relying on them.

### 7.3. Energy-histogram overlap between adjacent temperatures

**Script:** `scripts/analysis/check_energy_histogram_overlap.py`

This diagnostic reads the second comma-separated column (total energy) of each raw `mefile.txt`, converts it to energy per spin, and compares histograms at adjacent measured temperatures for each lattice size. A pair of histograms uses common bin edges. The reported overlap is

```text
hist_overlap = sum_i min(p_i, q_i)
```

where `p_i` and `q_i` are the normalized histogram-bin probabilities. The default is 25 bins, and the script plots the six adjacent-temperature pairs with the lowest overlap at each size.

For example:

```bash
python3 scripts/analysis/check_energy_histogram_overlap.py \
    --root results/square_cell_simulations \
    --output data/analysis_csv/diagnostics/square/cell/energy_overlap \
    --bins 25 --plots-per-size 6
```

The outputs include:

```text
adjacent_pairs.csv
summary_by_size.csv
overlap_vs_T_L<L>.png
size_<L>/worst_pair_<index>.png
```

`adjacent_pairs.csv` records each temperature pair, its overlap, and the separation between mean energies in units of the pooled standard deviation. `summary_by_size.csv` reports minimum and median overlaps and the fraction of adjacent pairs below fixed overlap thresholds. This is a **bin-dependent, finite-sample diagnostic**, not a proof that multihistogram reweighting is valid or invalid.

### 7.4. Cross-temperature correlations

**Script:** `scripts/analysis/check_cross_temperature_correlations_multi.py`

At each fixed lattice size, this script compares recorded series at different temperatures in the prepared `conf-<T>-<run>.dat` files. It computes Pearson correlations using the original saved-row ordering for energy, signed/absolute magnetization, signed/absolute Euler characteristic, and signed/absolute `D`.

It also evaluates the maximum absolute correlation across all temperature pairs and across adjacent temperature pairs. The null calibration independently circularly shifts the time series, preserving each individual series' internal temporal ordering up to the shift. The default is 500 null repetitions.

For example:

```bash
python3 scripts/analysis/check_cross_temperature_correlations_multi.py \
    --root data/eul2d/square/cell \
    --lattice square \
    --output data/analysis_csv/diagnostics/square/cell/cross_temperature \
    --expected-n 200 --null-repeats 500
```

The main outputs are:

```text
summary.csv
size_<L>/correlations_<observable>.csv
size_<L>/pairs_<observable>.csv
size_<L>/heatmap_<observable>.png
```

The script requires an unambiguous temperature-to-file mapping within each size; multiple runs at the same temperature are not automatically combined. Use `--expected-n 0` if measurement counts differ from the default 200 but remain consistent within each analyzed size. The resulting p-values are exploratory tests of **zero-lag linear correlation in saved-row order**; a nonsignificant result does not establish independence of the underlying random-number streams.

### 7.5. Cross-size correlations

**Script:** `scripts/analysis/cross_size_correlations.py`

This script performs a complementary check at **fixed temperature**, comparing time series from different lattice sizes. It uses Pearson correlations and circular-shift null calibrations, and provides results at each temperature, aggregated by observable across temperatures, and aggregated across all tested observables.

For example:

```bash
python3 scripts/analysis/cross_size_correlations.py \
    --root data/eul2d/square/cell \
    --lattice square \
    --output data/analysis_csv/diagnostics/square/cell/cross_size \
    --expected-n 200 --null-repeats 500
```

The principal outputs are:

```text
summary.csv
global_summary.csv
overall_summary.csv
T_<temperature>/correlations_<observable>.csv
T_<temperature>/pairs_<observable>.csv
T_<temperature>/heatmap_<observable>.png
```

By default, heatmaps are produced for every tenth temperature; use `--plot-step 1` to plot every temperature or `--plot-step 0` to disable heatmaps. If several files correspond to one `(L, T)` pair, use `--run` to select a run rather than combining them automatically. The default expected time-series length is 200, and `--expected-n 0` permits a common different length at a given temperature.

The cross-size and cross-temperature checks probe synchronized **linear fluctuations**, not correlations between sample means or the full joint distribution. The tests cannot establish statistical independence, including independence of the pseudorandom-number streams when seeds are reused.

### 7.6. Interpretation and relationship to the main analysis

These five scripts are supporting diagnostics, not additional steps required to redraw the publication figures from the processed CSVs.

- The main multihistogram stage produces bootstrap realizations used to locate susceptibility maxima.
- The leading FSS analysis uses the bootstrap standard deviation of pseudocritical peak locations in weighted fits.
- The square-lattice subleading FSS analysis uses paired resampling of those peak locations to quantify parameter uncertainty and contrasts.
- The raw temperature-dependent observable plots use nonoverlapping block-jackknife errors, with their own autocorrelation-based block-size selection.
- The simulation runners use `Ndecorr = 10` between saved measurements. In the current implementation, the seed file is reused rather than updated with an independent seed for every run. Neither a lack of detected cross-stream correlation nor reuse of the same seed should be described as evidence of independent random-number initialization.

The diagnostics require raw or prepared simulation data, whereas the processed data included in the repository are sufficient for the documented figure-reproduction commands. The example output directories under `data/analysis_csv/diagnostics/` are choices for organizing newly generated diagnostics; they are **not required input paths** for the principal plotting scripts.
