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

## Finite-size-scaling figure

The finite-size-scaling figure compares pseudocritical temperatures obtained
from the Euler-characteristic and magnetization observables.

The plotting script is:

    scripts/plotting/plot_fss_tc_two_observables.py

It reads the previously generated pseudocritical-temperature summaries:

    data/analysis_csv/<lattice_type>/pseudocritical_ec.csv
    data/analysis_csv/<lattice_type>/pseudocritical_m.csv

The fitted finite-size-scaling model is:

    Tc(L) = Tc_inf + A * L^(-1/nu)

The EC- and magnetization-based pseudocritical temperatures are fitted
independently.

The default lower lattice-size limits used in the fits are:

    square:
        EC: L >= 96
        M:  L >= 96

    triangular:
        EC: L >= 96
        M:  L >= 64

These limits are analysis choices and can be overridden explicitly with:

    --lmin-ec
    --lmin-m

For example:

    python scripts/plotting/plot_fss_tc_two_observables.py \
      --lattice-type square \
      --lmin-ec 128 \
      --lmin-m 128

By default, the exponent nu is fitted independently for each observable.
It can optionally be fixed using:

    --fix-nu-ec
    --fix-nu-m

The horizontal plotting coordinate is:

    x = L^(-1/nu_ref)

with nu_ref = 1 by default, corresponding to x = 1/L.

The default figure commands are:

    python scripts/plotting/plot_fss_tc_two_observables.py \
      --lattice-type square

and:

    python scripts/plotting/plot_fss_tc_two_observables.py \
      --lattice-type triangular

The generated figures are written to:

    results/figures/square_fss_Tc_EC_vs_M.png
    results/figures/triangular_fss_Tc_EC_vs_M.png

The figure uses the existing pseudocritical-temperature CSV files directly;
no additional figure-level data file is required.

## FSS stability with fixed maximum lattice size

The stability of the finite-size-scaling estimates is examined by varying the
minimum lattice size included in the fit while keeping the maximum lattice size
fixed.

The analysis script is:

    scripts/analysis/compute_fss_fixed_lmax.py

It reads the previously generated pseudocritical-temperature summaries:

    data/analysis_csv/<lattice_type>/pseudocritical_ec.csv
    data/analysis_csv/<lattice_type>/pseudocritical_m.csv

For each observable, the model

    Tc(L) = Tc_inf + A * L^(-1/nu)

is fitted repeatedly while scanning L_min.

Two fixed upper limits are considered:

    L_max = 3072
    L_max = 2048

corresponding to the largest and second-largest lattice sizes in the dataset.

Fits are performed only when at least six lattice sizes remain in the selected
range. AICc and additional fit diagnostics are also recorded, although the
stability figures themselves use the fitted values of Tc_inf and nu together
with their uncertainties.

The resulting tables are stored in:

    data/table_csv/<lattice_type>/fixed_lmax/

with files:

    lmin_ec_lmax_max.csv
    lmin_m_lmax_max.csv
    lmin_ec_lmax_second.csv
    lmin_m_lmax_second.csv

The plotting script is:

    scripts/plotting/plot_fss_lmin_stability.py

It generates plots of fitted Tc_inf and nu as functions of L_min for both
Euler-characteristic and magnetization observables.

For the square lattice, the theoretical reference values are:

    nu = 1
    Tc = 2 / ln(1 + sqrt(2))

For the triangular lattice, the theoretical reference values are:

    nu = 1
    Tc = 4 / ln(3)

The plotting script processes square and triangular datasets independently. If
one lattice dataset is absent, that lattice is skipped without stopping the
generation of figures for the other lattice.

The square-lattice output files corresponding to the publication panels are:

    results/figures/Tc_vs_Lmin_PBC_Lmax_max.png
    results/figures/Tc_vs_Lmin_PBC_Lmax_second.png
    results/figures/nu_vs_Lmin_PBC_Lmax_max.png
    results/figures/nu_vs_Lmin_PBC_Lmax_second.png

The corresponding triangular-lattice figures are generated with separate
triangular output names.

## Distribution-shape and effective-resolution analysis

The final analysis compares the distributions of magnetization and the
symmetrized Euler characteristic at a selected temperature close to the
critical region.

The analysis script is:

    scripts/analysis/compute_distribution_resolution.py

It reads the raw simulation files:

    mefile.txt
    epfile.txt
    enfile.txt

for each lattice size at a user-selected temperature.

The symmetrized Euler characteristic is defined as:

    EC_sym = (EN - EP) / 2

For the distribution-shape comparison, magnetization and EC_sym are
standardized independently using their sample mean and sample standard
deviation. A two-sample Kolmogorov-Smirnov statistic is then computed between
the standardized distributions.

The resulting table is stored in:

    data/analysis_csv/<lattice_type>/distribution_resolution/
        distribution_shape_tests.csv

It contains, among other quantities:

    L
    n_samples
    M_mean
    M_std
    M_skewness
    M_excess_kurtosis
    EC_mean
    EC_std
    EC_skewness
    EC_excess_kurtosis
    KS_statistic_standardized
    KS_pvalue_standardized

The same raw observables are also used to quantify their effective discrete
resolution.

For an empirical distribution with probabilities p_i, the Shannon entropy is:

    H = -sum_i p_i log(p_i)

and the effective number of observed values is:

    K_eff = exp(H)

The corresponding table is stored in:

    data/analysis_csv/<lattice_type>/distribution_resolution/
        entropy_effective_resolution.csv

It contains the number of unique observed values, Shannon entropy, effective
number of values, and EC-to-magnetization resolution ratios.

For the square-lattice publication figures, the analysis is evaluated at:

    T = 2.26900

For example:

    python scripts/analysis/compute_distribution_resolution.py \
      --results-dir <square_simulation_results_directory> \
      --lattice-type square \
      --temperature 2.26900

The temperature is supplied explicitly through the command line and is not
hard-coded into the analysis.

### Distribution-resolution figures

The plotting script is:

    scripts/plotting/plot_distribution_resolution.py

It reads only the processed CSV files and generates:

    results/figures/KS_distance_standardized_vs_L.png
    results/figures/resolution_ratios_EC_over_M_vs_L.png

The first figure shows the Kolmogorov-Smirnov distance between the standardized
magnetization and EC_sym distributions as a function of lattice size.

The second figure compares:

    K_EC / K_M

and:

    K_eff,EC / K_eff,M

as functions of lattice size. A horizontal reference line at one indicates
equal effective resolution for the two observables.

The complete workflow for this analysis is:

    raw simulation outputs
            |
            v
    scripts/analysis/compute_distribution_resolution.py
            |
            v
    distribution_shape_tests.csv
    entropy_effective_resolution.csv
            |
            v
    scripts/plotting/plot_distribution_resolution.py
            |
            v
    KS_distance_standardized_vs_L.png
    resolution_ratios_EC_over_M_vs_L.png
