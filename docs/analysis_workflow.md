# Analysis and figure reproduction

This document describes the analysis steps used to transform simulation outputs
into lightweight figure-level data and publication figures.

The analysis is separated from the simulation and multihistogram workflows.
Raw simulation outputs are first reduced to the quantities required for a
specific analysis. These processed data are stored as CSV files, which are then
used directly by the plotting scripts.

This separation makes it possible to reproduce publication figures without
rerunning the full simulations.

## Temperature-dependent observables

The temperature-dependent magnetization, energy, and Euler-characteristic
observables are computed directly from the raw square- and triangular-lattice
simulation outputs.

The analysis script is:

    scripts/analysis/compute_observables_vs_temperature.py

It reads the simulation files:

    mefile.txt
    epfile.txt
    enfile.txt

for a selected lattice size and computes:

- <|M|>
- <E>
- <|EC_sym|>, where EC_sym = (chi_neg - chi_pos) / 2
- <EC_avg>, where EC_avg = (chi_neg + chi_pos) / 2

Uncertainties are estimated using the jackknife standard error of the mean.

The default lattice size used for these figures is:

    L = 3072

and can be changed with the `--size` argument.

## Input data

By default, the analysis script expects simulation outputs in:

    results/square_simulations/
    results/triangular_simulations/

with the directory structure:

    results/<lattice_type>/T_<temperature>/size_<L>/

Alternative locations can be supplied explicitly:

    python scripts/analysis/compute_observables_vs_temperature.py \
      --square-results <square_results_directory> \
      --triangular-results <triangular_results_directory>

This allows the analysis to be performed using separately distributed
simulation data without copying the full dataset into the repository.

## Figure-level data

The analysis produces lightweight CSV files in:

    data/figure_csv/

The current outputs are:

    square_mag_vs_T.csv
    square_energy_vs_T.csv
    square_euler_sym_vs_T.csv
    square_ec_avg_vs_T.csv

    triangular_mag_vs_T.csv
    triangular_energy_vs_T.csv
    triangular_euler_sym_vs_T.csv
    triangular_ec_avg_vs_T.csv

Each file contains the columns:

    lattice_type
    L
    T
    observable
    mean
    jackknife_se
    n_samples

These CSV files contain the processed numerical data required to reproduce the
corresponding figures.

## Generating the figures

Figures are generated from the figure-level CSV files with:

    scripts/plotting/plot_observables_vs_temperature.py

Run:

    python scripts/plotting/plot_observables_vs_temperature.py

to generate the figures in:

    results/figures/

The generated files are:

    2d_mag_vs_T_abs.png
    2d_energy_vs_T.png
    2d_euler_sym_vs_T_abs.png
    2d_ec_avg_vs_T.png

    tri_mag_vs_T_abs.png
    tri_energy_vs_T.png
    tri_euler_sym_vs_T_abs.png
    tri_ec_avg_vs_T.png

The plotting script reads only the processed CSV files. It does not recompute
observables or statistical uncertainties.

## Workflow summary

The current analysis pipeline is:

    raw simulation outputs
            |
            v
    scripts/analysis/compute_observables_vs_temperature.py
            |
            v
    data/figure_csv/*.csv
            |
            v
    scripts/plotting/plot_observables_vs_temperature.py
            |
            v
    results/figures/*.png

Additional publication figures and their corresponding analysis steps will be
documented here as they are added to the reproducibility workflow.

## Pseudocritical temperatures and finite-size scaling

Pseudocritical temperatures are extracted from the multihistogram bootstrap
outputs stored as:

    results_L.dat

The extraction script is:

    scripts/analysis/compute_pseudocritical_temperatures.py

For each lattice size and bootstrap realization, the pseudocritical temperature
is identified from the maximum of the selected susceptibility-like observable.

Two variants are analyzed:

- `ec`: Euler-characteristic-based observable
- `m`: magnetization-based observable

The extracted bootstrap-level pseudocritical temperatures are written to:

    data/analysis_csv/<lattice_type>/pseudocritical_bootstrap_ec.csv
    data/analysis_csv/<lattice_type>/pseudocritical_bootstrap_m.csv

The corresponding size-dependent summaries are written to:

    data/analysis_csv/<lattice_type>/pseudocritical_ec.csv
    data/analysis_csv/<lattice_type>/pseudocritical_m.csv

The summary files contain:

    lattice_type
    variant
    L
    Tc
    Tc_std
    N_boot

where `Tc` is the mean pseudocritical temperature across bootstrap realizations
and `Tc_std` is its sample standard deviation.

### Finite-size-scaling fits

Finite-size scaling is performed using:

    scripts/analysis/compute_fss_lmin_aic.py

The script reads a pseudocritical-temperature summary CSV and fits

    Tc(L) = Tc_inf + A * L^(-1/nu)

for a sequence of minimum lattice sizes `L_min`.

For each fit, the output includes:

    L_min
    n_points
    Tc_inf
    dTc_inf
    A
    dA
    nu
    dnu
    chi2
    dof
    chi2_red
    p_value
    AIC
    dnu_sigma
    dTc_inf_sigma
    delta_AIC

A finite-sample AIC correction (AICc) is used by default.

The resulting tables are stored separately for square and triangular lattices,
for example:

    data/table_csv/square/aic_lmin_ec.csv
    data/table_csv/square/aic_lmin_m.csv

    data/table_csv/triangular/aic_lmin_ec.csv
    data/table_csv/triangular/aic_lmin_m.csv

The corresponding analysis pipeline is:

    multihistogram results_L.dat
            |
            v
    scripts/analysis/compute_pseudocritical_temperatures.py
            |
            v
    data/analysis_csv/<lattice_type>/pseudocritical_*.csv
            |
            v
    scripts/analysis/compute_fss_lmin_aic.py
            |
            v
    data/table_csv/<lattice_type>/aic_lmin_*.csv
