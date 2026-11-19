
# Ising reproducibility

This repository contains the simulation, multihistogram-reweighting, statistical-analysis, and plotting code associated with the study of topological observables in the two-dimensional Ising model.

## Associated publication

**Topological signatures of criticality: the Ising case**

Tak-Shing Chan, Pawe³ D³otko, Jacek Gulgowski, Biagio Lucini, Mateusz Mas³owski, and Davide Vadacchino.

Preprint: [arXiv:2611.23586](https://arxiv.org/abs/2611.23586) (submitted 17 November 2026).

Journal status: submitted to *Physical Review E*.

The published-analysis code version and associated data release will be linked here once permanent archive identifiers are available.

## Scope

The code supports square and triangular lattices and two representations of spin domains (`cell` and `vertex`). The analysis used for the associated paper focuses on:

- **Square lattice:** `square/cell`
- **Triangular lattice:** `triangular/vertex`

The principal observables include magnetization, energy, the symmetrized Euler characteristic, and a related observable `d`. The analysis uses multihistogram reweighting to estimate pseudocritical temperatures, performs finite-size-scaling (FSS) fits, and compares subleading corrections on the square lattice.

The C/C++ simulation and reweighting code can be rebuilt locally. Large raw simulation and multihistogram datasets are **not** included in the Git repository. Selected processed CSVs and figures **are** tracked, so many plots can be regenerated without rerunning production simulations.

## Repository layout

```text
src/
  simulation/             C simulation code and Makefile
  multihistogram/         C++ reweighting code and Makefile

scripts/
  simulation/             Square/triangular simulation runners
  multihistogram/         Temperature grids, input preparation, reweighting
  analysis/               Statistical analysis and FSS scripts
  plotting/               Plot-generation scripts

docs/
  simulation_workflow.md
  multihistogram_workflow.md
  analysis_workflow.md
  software_requirements.md

data/
  analysis_csv/           Processed observables, pseudocritical peaks,
                          resolution data and paired FSS results
  table_csv/              FSS scan tables
  analysis_parameters/    Parameters used in the analysis
  temperatures/           Generated temperature-grid files (not tracked)
  eul2d/                 Generated multihistogram input (not tracked)

results/
  figures/                Selected generated figures (tracked)
  multihistogram/         Generated reweighting results (not tracked)

requirements.txt          Python dependencies
```

Generated directories are populated by the relevant scripts. The repository already includes some processed data products; it is not necessary to rerun the simulation to regenerate figures from these files.

## Requirements

A Unix-like environment with Bash, GNU Make, GCC/G++, `awk`, and standard shell tools. Python dependencies are listed in `requirements.txt` (NumPy, pandas, Matplotlib, and SciPy).

From the repository root:

```bash
python3 -m pip install -r requirements.txt
make -C src/simulation
make -C src/multihistogram
```

The simulation Makefile generates executables for the lattice sizes `64 96 128 192 256 384 512 768 1024 1536 2048 3072`.

Its executable names have the form `ising2d-cell-L`, `ising2d-vertex-L`, `isingtr-cell-L`, and `isingtr-vertex-L`.

The multihistogram executable is `src/multihistogram/Rw`.

See [software requirements](docs/software_requirements.md) for build details.

## Reproduce figures using the included CSV files

This is the shortest route when only the existing processed-data plots are needed. Run from the repository root:

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

Figures are written under `results/figures/<lattice>/<representation>/`. Temperature-dependent observable plots are stored under the additional subdirectory `observables_vs_temperature/size_3072/`.

The plotting scripts consume processed CSVs from `data/analysis_csv/` and `data/table_csv/`; they **do not** repeat the simulations or multihistogram calculations. The plots of temperature-dependent observables use block-jackknife standard errors saved in the corresponding CSVs.

## Full workflow, from simulations to multihistogram outputs

The following examples use the two lattice/representation combinations analyzed in the paper. They may require substantial computing resources; the simulation runners can be dispatched through a cluster scheduler.

### 1. Run simulations

Build the simulation executables first, as shown above. Then:

```bash
bash scripts/simulation/run_square_2d_simulations.sh \
    src/simulation results/square_cell_simulations cell

bash scripts/simulation/run_triangular_2d_simulations.sh \
    src/simulation results/triangular_vertex_simulations vertex
```

Current runner defaults are:

| Lattice | Start temperature | End temperature | Step | Configurations |
|---|---:|---:|---:|---:|
| Square | 2.22 | 2.32 | 0.0005 | 200 |
| Triangular | 3.59 | 3.72 | 0.0005 | 200 |

The optional positional arguments after `cell`/`vertex` override the start temperature, end temperature, step, and number of configurations, in that order. Each simulation results directory uses a `seed.txt` file, initialized to `123456` if absent and subsequently reused.

### 2. Generate reweighting temperature grids

```bash
bash scripts/multihistogram/generate_square_reweighting_temperatures.sh
bash scripts/multihistogram/generate_triangular_reweighting_temperatures.sh
```

These generate `data/temperatures/square/temperatures_L.txt` and `data/temperatures/triangular/temperatures_L.txt`, with temperature windows depending on lattice size.

### 3. Prepare multihistogram inputs

The input-preparation runner accepts four positional arguments:

`<lattice_type> <representation> <temperatures_dir> <simulation_results_dir>`

It creates the corresponding `data/eul2d/<lattice>/<representation>/` directory.

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
    square cell \
    data/temperatures/square \
    results/square_cell_simulations

bash scripts/multihistogram/prepare_eul2d_input.sh \
    triangular vertex \
    data/temperatures/triangular \
    results/triangular_vertex_simulations
```

### 4. Run multihistogram reweighting

Build `src/multihistogram/Rw` first.

The runner accepts:

`<lattice_type> <representation> <mh_executable_dir> [refinement_factor]`

The default refinement factor is `50`.

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
    square cell src/multihistogram

bash scripts/multihistogram/run_multihistogram_all.sh \
    triangular vertex src/multihistogram
```

The reweighting results are written to:

```text
results/multihistogram/square/cell/results_L.dat
results/multihistogram/triangular/vertex/results_L.dat
```

Other lattice/representation combinations can be run using the corresponding `cell`/`vertex` arguments and simulation outputs; they are not the focus of the paper's final analysis.

## Statistical analysis

The main analysis stages are described below. For each stage, use the corresponding raw or reweighted data from the desired lattice/representation combination.

### Pseudocritical temperatures

For each bootstrap realization and lattice size, `compute_pseudocritical_temperatures.py` extracts the location of the susceptibility peak for `ec`, `m`, and `d`. It produces paired bootstrap samples and per-size summary CSVs.

```bash
for combo in "square cell" "triangular vertex"; do
    read -r LAT REP <<< "$combo"
    python3 scripts/analysis/compute_pseudocritical_temperatures.py \
        --input-dir "results/multihistogram/$LAT/$REP" \
        --lattice-type "$LAT" \
        --representation "$REP" \
        --output-dir data/analysis_csv
done
```

Each combination produces `pseudocritical_{ec,m,d}.csv` and `pseudocritical_bootstrap_{ec,m,d}.csv` under `data/analysis_csv/<lattice>/<representation>/`.

### Leading FSS and stability scans

The leading fit uses

$$
T_c(L) = T_c(\infty) + A L^{-1/\nu}.
$$

Generate FSS scans over the minimum lattice size:

```bash
for combo in "square cell" "triangular vertex"; do
    read -r LAT REP <<< "$combo"
    for OBS in ec m d; do
        python3 scripts/analysis/compute_fss_lmin.py \
            --input "data/analysis_csv/$LAT/$REP/pseudocritical_${OBS}.csv" \
            --output "data/table_csv/$LAT/$REP/lmin_${OBS}.csv"
    done

    python3 scripts/analysis/compute_fss_fixed_lmax.py \
        --lattice-type "$LAT" \
        --representation "$REP" \
        --input-dir data/analysis_csv \
        --output-dir data/table_csv \
        --min-points 6
done
```

The fixed-`L_max` scan uses the largest and second-largest lattice sizes available and supplies the CSVs used by `plot_fss_lmin_stability.py`. Fit quality is evaluated using chi-square goodness-of-fit probabilities.

### Subleading square-lattice corrections

For the `square/cell` dataset, the correction fit fixes `nu = 1` and uses

$$
T_c(L) = T_c(\infty) + \frac{a}{L} + \frac{b}{L^2}.
$$

```bash
python3 scripts/analysis/paired_fss_corrections.py \
    --lattice-type square \
    --representation cell \
    --input-base data/analysis_csv \
    --output-base data/analysis_csv/fss_corrections \
    --lmin 64 --lmax 3072 \
    --n-resamples 10000 --seed 123456 --ci 95
```

This script fits the mean pseudocritical temperatures using their standard deviations as weights. Paired bootstrap samples within each lattice size are then used to estimate parameter uncertainties and differences in the subleading coefficient `b` between the observables.

The outputs include:

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

### Temperature-dependent observables

The following example calculates raw-simulation observables at `L=3072`, using automatic block-size selection for the block jackknife, then creates the four associated plots:

```bash
python3 scripts/analysis/compute_observables_vs_temperature.py \
    --results-dir results/square_cell_simulations \
    --lattice-type square --representation cell --size 3072 \
    --output-dir data/analysis_csv \
    --tau-multiplier 10 --min-blocks 10

python3 scripts/plotting/plot_observables_vs_temperature.py \
    --lattice-type square --representation cell --size 3072
```

The data are stored under `data/analysis_csv/square/cell/observables_vs_temperature/size_3072/`.

Replace the lattice, representation and results directory to analyze the triangular-lattice data. The `jackknife_se` column contains block-jackknife errors; the script does not use the paired FSS bootstrap for these raw-observable plots.

### Distribution and effective resolution

The distribution-resolution workflow reads raw spin-observable time series at an explicitly selected simulation temperature and produces `distribution_shape.csv` and `entropy_effective_resolution.csv`.

The current figure plots the relative unique-state and entropy-based effective-resolution measures, without Kolmogorov-Smirnov plots.

For example, if the simulation directory `T_2.26900` exists:

```bash
python3 scripts/analysis/compute_distribution_resolution.py \
    --results-dir results/square_cell_simulations \
    --lattice-type square --representation cell \
    --temperature 2.26900 --output-dir data/analysis_csv

python3 scripts/plotting/plot_distribution_resolution.py \
    --lattice-type square --representation cell
```

Additional scripts under `scripts/analysis/` provide autocorrelation, cross-temperature/cross-size correlation, and histogram-overlap diagnostics.

## Data availability and reproducibility limits

The Git repository contains source code, scripts, selected analysis CSVs, selected FSS table CSVs, and figure PNGs.

It does **not** contain the complete raw time series, prepared multihistogram input, or complete reweighting outputs needed to rerun every analysis step from scratch.

**External data archive:** link to be added after archival release.

For full reproduction, place the external simulation/reweighting data in the directory layouts described above (or supply matching directories to scripts that accept path arguments).

For figure reproduction, the processed CSV files tracked in this repository are sufficient for the supported plotting scripts.

Production simulations may be expensive at the largest lattice sizes and can be scheduled or parallelized according to the local computing environment.

For the exact software versions, see `requirements.txt` and the compiler/build documentation.

## Documentation and citation

- [Simulation workflow](docs/simulation_workflow.md)
- [Multihistogram workflow](docs/multihistogram_workflow.md)
- [Analysis workflow](docs/analysis_workflow.md)
- [Software requirements](docs/software_requirements.md)

When citing the scientific results, please cite [the associated preprint](https://arxiv.org/abs/2611.23586), or the journal article once available.

A `CITATION.cff` file and repository/data archive identifiers can be added as part of the final release.

**Licensing:** The repository does not currently contain a top-level `LICENSE` file. Licensing terms should be clarified before the archived release.
