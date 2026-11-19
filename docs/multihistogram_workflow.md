# Multihistogram workflow

## Overview

This document describes how to prepare simulation measurements and run multihistogram reweighting for the two-dimensional Ising model. The workflow supports square and triangular lattices, each with `cell` and `vertex` representations. The lattice/representation combinations analyzed in the associated publication are:

- `square/cell`
- `triangular/vertex`

The workflow has four stages:

1. Generate lattice-size-dependent temperature files from the tracked parameter table.
2. Convert raw simulation measurements into the ten-column multihistogram input format.
3. Run the C++ reweighting executable with circular moving block bootstrap resampling.
4. Use the reweighted bootstrap results in the pseudocritical-temperature and finite-size-scaling analyses.

Unless otherwise specified, run the commands below from the repository root.

## Build the multihistogram executable

```bash
make -C src/multihistogram
```

This builds:

```text
src/multihistogram/Rw
```

The Makefile uses `g++` and the C++ sources in `src/multihistogram/`, including the local autocorrelation and RANLUX implementations. See [Software requirements](software_requirements.md) for additional details.

## 1. Generate temperature files

The temperature windows are recorded in:

```text
data/analysis_parameters/reweighting_windows.csv
```

Each row specifies `lattice_type`, `L`, `T_min`, `T_max`, and `T_step`. Windows vary with lattice size. The current table uses `T_step = 0.0005` throughout; the step is read from the CSV rather than hard-coded into the generator.

Generate temperature files for both lattices:

```bash
bash scripts/multihistogram/generate_square_reweighting_temperatures.sh
bash scripts/multihistogram/generate_triangular_reweighting_temperatures.sh
```

The generated paths are:

```text
data/temperatures/square/temperatures_L.txt
data/temperatures/triangular/temperatures_L.txt
```

Each temperature file contains two whitespace-separated columns:

```text
T IDX
```

For example, the square-lattice file for `L=128` begins with:

```text
2.2650 0
2.2655 1
2.2660 2
```

Here `IDX` is the integer identifier used in the corresponding input filename. The generators format temperatures to four decimal places.

Both generation scripts accept optional positional arguments:

```text
[parameter_csv] [output_directory]
```

For example:

```bash
bash scripts/multihistogram/generate_square_reweighting_temperatures.sh \
    data/analysis_parameters/reweighting_windows.csv \
    data/temperatures/square
```

The temperature windows are analysis inputs specified by the CSV. They are not selected automatically by the generator.

## 2. Prepare multihistogram input files

The input-preparation script is:

```text
scripts/multihistogram/prepare_eul2d_input.sh
```

It requires four positional arguments:

```text
<lattice_type> <representation> <temperatures_dir> <simulation_results_dir>
```

For the square-lattice cell representation:

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
    square cell \
    data/temperatures/square \
    results/square_cell_simulations
```

For the triangular-lattice vertex representation:

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
    triangular vertex \
    data/temperatures/triangular \
    results/triangular_vertex_simulations
```

The same script also supports `square/vertex` and `triangular/cell`, provided corresponding simulation outputs exist.

The preparation script reads raw measurements from:

```text
<simulation_results_dir>/T_<temperature>/size_<L>/
    mefile.txt
    epfile.txt
    enfile.txt
    fpfile.txt
    fnfile.txt
```

Simulation directories usually format temperatures to five decimal places. The preparation script matches them against four-decimal temperatures from `temperatures_L.txt`.

For each available temperature and lattice size, it writes:

```text
data/eul2d/<lattice_type>/<representation>/size_<L>/conf-<T>-<IDX>.dat
```

For example:

```text
data/eul2d/square/cell/size_128/conf-2.2650-0.dat
data/eul2d/triangular/vertex/size_64/conf-3.6800-0.dat
```

### Input columns

Each `conf-<T>-<IDX>.dat` file contains ten whitespace-separated columns, with one row per recorded configuration:

| Column | Quantity | Source |
|---|---|---|
| 1 | Magnetization `M` | First field of `mefile.txt` |
| 2 | Energy `E` | Second field of `mefile.txt` |
| 3 | Positive-spin Euler characteristic `chi_plus` | `epfile.txt` |
| 4 | Negative-spin Euler characteristic `chi_minus` | `enfile.txt` |
| 5-7 | Positive-spin cell counts | First three fields of `fpfile.txt` |
| 8-10 | Negative-spin cell counts | First three fields of `fnfile.txt` |

The files `bpfile.txt` and `bnfile.txt` are not used by this preparation step.

The preparation script checks whether its five required input files exist and are nonempty. If their lengths differ, it uses the smallest row count and therefore truncates longer files. Missing temperature directories or measurements are skipped and reported. Review the output before reweighting: the C++ executable expects a prepared file for every temperature listed in the relevant temperature file.

## 3. Run multihistogram reweighting

The runner is:

```text
scripts/multihistogram/run_multihistogram_all.sh
```

It takes three required arguments and one optional argument:

```text
<lattice_type> <representation> <mh_executable_dir> [refinement_factor]
```

For the combinations analyzed in the publication:

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
    square cell src/multihistogram

bash scripts/multihistogram/run_multihistogram_all.sh \
    triangular vertex src/multihistogram
```

The runner determines input and output locations from the lattice and representation. It also passes `lattice_type` (`square` or `triangular`) to the `Rw` executable as its final command-line argument:

```text
Rw temperatures_file L dimension results_file T_min T_max n_points output_file input_data_dir lattice_type
```

Both `input_data_dir` and `lattice_type` are required by the updated executable. The runner supplies these arguments automatically; users invoking `Rw` directly must supply them explicitly.

```text
Temperature files:
    data/temperatures/<lattice_type>/temperatures_L.txt

Prepared input:
    data/eul2d/<lattice_type>/<representation>/size_L/

Output:
    results/multihistogram/<lattice_type>/<representation>/
```

It processes the supported sizes:

```text
64 96 128 192 256 384 512 768 1024 1536 2048 3072
```

A size is skipped if its temperature file or prepared input directory is missing. Once reweighting begins for a size, a missing or malformed individual configuration file causes the C++ executable to fail.

### Temperature-grid refinement

For each lattice size, the runner determines the number of reweighting output points as:

```text
POINTS = (N_temperatures - 1) * refinement_factor + 1
```

The default `refinement_factor` is `50`; it can be supplied as the fourth argument:

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
    square cell src/multihistogram 100
```

**Important:** the C++ executable distributes reweighting points uniformly in inverse temperature, `beta = 1/T`, between `1/T_max` and `1/T_min`. Consequently, the output temperature spacing is not exactly uniform. For input spacings of `0.0005` and refinement factor `50`, the typical output temperature spacing is approximately `0.00001`.

### Circular moving block bootstrap

The current `src/multihistogram/main.cpp` uses circular moving block bootstrap resampling on the recorded configurations at each simulation temperature.

For every input time series, integrated autocorrelation times are estimated using `AutoCorr` for seven observables:

- Energy `E`
- Absolute magnetization `|M|`
- Absolute symmetrized Euler characteristic `|EC_sym|`
- Absolute derived observable `|d|`
- Squares of the last three absolute observables

The maximum estimated autocorrelation time determines a common block length for energy and all observables at that temperature:

```text
block_length = ceil(block_factor * max(0.5, tau_max))
```

The default `block_factor` is `10`. Blocks start at uniformly sampled positions, wrap around the end of the recorded series, and are sampled with replacement until one bootstrap sample of the original length has been assembled. The final block is truncated if necessary.

The same resampled configuration indices are used for energy and all observables at a given temperature, preserving their within-configuration relationships. For a given bootstrap replicate, blocks are sampled separately for each temperature. The procedure does not itself establish statistical independence of the original simulation streams across temperatures or lattice sizes.

The default number of bootstrap replicates is `200`. The executable accepts three environment variables:

| Variable | Default | Purpose |
|---|---|---|
| `ISING_BLOCK_FACTOR` | `10` | Factor multiplying the maximum estimated autocorrelation time when determining block length |
| `ISING_BOOTSTRAP_REPS` | `200` | Number of bootstrap replicates |
| `ISING_BOOTSTRAP_SEED` | `123456` | Seed for the bootstrap random-number generator |

For an explicit run with these defaults:

```bash
ISING_BLOCK_FACTOR=10 \
ISING_BOOTSTRAP_REPS=200 \
ISING_BOOTSTRAP_SEED=123456 \
bash scripts/multihistogram/run_multihistogram_all.sh \
    square cell src/multihistogram
```

The bootstrap seed is separate from the seed file used by the Monte Carlo simulation programs. The `AutoCorr` estimates are expressed in intervals between **saved measurements**, not individual Metropolis or Wolff updates. The bootstrap executable does not enforce a minimum number of nonoverlapping blocks; separate diagnostic scripts can assess the feasibility and sensitivity of block-size choices.

### Lattice-dependent definition of `d`

The current shared executable takes `lattice_type` as a required command-line argument and uses

```text
EC_sym = (chi_plus - chi_minus) / 2

d = M + 2 * EC_sym    (square)
d = M + EC_sym        (triangular)
```

In `main.cpp`, the coefficient is selected from the supplied lattice type before the input measurements are processed. The observable passed to reweighting is `|d| / L^2`. This agrees with the lattice-dependent formula used in the autocorrelation and cross-stream diagnostic scripts. The symmetrized Euler characteristic and magnetization observables are treated in the same way for both lattice types.

## 4. Output files and post-processing

For each lattice size, the runner creates two files:

```text
results/multihistogram/<lattice_type>/<representation>/results_L.dat
results/multihistogram/<lattice_type>/<representation>/output2_L.dat
```

Examples:

```text
results/multihistogram/square/cell/results_128.dat
results/multihistogram/square/cell/output2_128.dat

results/multihistogram/triangular/vertex/results_128.dat
results/multihistogram/triangular/vertex/output2_128.dat
```

Both files contain whitespace-separated numeric data without a header row. Each bootstrap replicate is identified by its zero-based `rep` index.

### Reweighted results: `results_L.dat`

There are nine columns:

| Column | Contents |
|---|---|
| 1 | Bootstrap replicate ID `rep` |
| 2 | Reweighted temperature `T` |
| 3 | Reweighted free-energy-related quantity `freeNRG` |
| 4 | Mean of `|EC_sym| / L^2` |
| 5 | `L^2` times its variance |
| 6 | Mean of `|M| / L^2` |
| 7 | `L^2` times its variance |
| 8 | Mean of `|d| / L^2` |
| 9 | `L^2` times its variance |

In this file, the variance columns do **not** include the additional factor `beta = 1/T`. The downstream pseudocritical-temperature script divides these columns by `T` when locating the peaks of the corresponding susceptibility.

### Estimates at simulated temperatures: `output2_L.dat`

There are eight columns:

| Column | Contents |
|---|---|
| 1 | Bootstrap replicate ID `rep` |
| 2 | Original simulation temperature `T` |
| 3 | Bootstrap mean of `|EC_sym| / L^2` |
| 4 | `L^2` times its sample variance (using the bootstrap sample moments) |
| 5 | Bootstrap mean of `|M| / L^2` |
| 6 | `L^2` times its sample variance (using the bootstrap sample moments) |
| 7 | Bootstrap mean of `|d| / L^2` |
| 8 | `L^2` times its sample variance (using the bootstrap sample moments) |

This file contains measurements evaluated at the input simulation temperatures, rather than on the refined reweighting grid.

### Extract pseudocritical temperatures

Once `results_L.dat` files are available, extract the bootstrap maxima for the `ec`, `m`, and `d` variants:

```bash
python3 scripts/analysis/compute_pseudocritical_temperatures.py \
    --input-dir results/multihistogram/square/cell \
    --lattice-type square \
    --representation cell \
    --output-dir data/analysis_csv

python3 scripts/analysis/compute_pseudocritical_temperatures.py \
    --input-dir results/multihistogram/triangular/vertex \
    --lattice-type triangular \
    --representation vertex \
    --output-dir data/analysis_csv
```

The processed files are written under `data/analysis_csv/<lattice_type>/<representation>/`. Subsequent FSS and figure-generation steps are documented in [Analysis and figure reproduction](analysis_workflow.md).

## Reproducibility and data availability

The parameter table, source code, build instructions, and workflow scripts are versioned in the repository. Temperature lists, prepared configuration files, and multihistogram outputs are generated by the workflow.

The complete production time series and multihistogram datasets are not stored in the Git repository. Full reproduction from raw measurements requires access to those datasets or rerunning the simulations.

The current C++ implementation reports autocorrelation estimates and selected block lengths during execution. Their reliability, and the sensitivity of uncertainties to block length, should be assessed using the statistical diagnostics described in [Analysis and figure reproduction](analysis_workflow.md).
