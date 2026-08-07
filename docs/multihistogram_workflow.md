# Multihistogram workflow

## Overview

The multihistogram analysis uses configuration files prepared separately for the square and triangular lattices.

The input data are organized as:

```text
data/eul2d/square/size_L/
data/eul2d/triangular/size_L/
```

where `L` is the linear lattice size of the `L x L` grid.

Each `size_L` directory contains files named:

```text
conf-T-IDX.dat
```

where:

* `T` is the simulation temperature formatted with four decimal places
* `IDX` is the integer index assigned to that temperature

Examples:

```text
data/eul2d/square/size_64/conf-2.2200-0.dat
data/eul2d/triangular/size_64/conf-3.6800-0.dat
```

## Building the multihistogram executable

From the repository root:

```bash
cd src/multihistogram
make
cd ../..
```

This builds the multihistogram executable:

```text
src/multihistogram/Rw
```

## Temperature files

Temperature files are stored separately for the two lattice types:

```text
data/temperatures/square/
data/temperatures/triangular/
```

The files are named:

```text
temperatures_L.txt
```

Each line has the format:

```text
T IDX
```

For example:

```text
2.2600 0
2.2605 1
2.2610 2
```

The reweighting temperature windows were selected manually around the relevant peak region for each lattice size.

### Square lattice

Square-lattice temperature files are generated with:

```bash
bash scripts/multihistogram/generate_square_reweighting_temperatures.sh
```

The default output directory is:

```text
data/temperatures/square/
```

### Triangular lattice

Triangular-lattice temperature files are generated with:

```bash
bash scripts/multihistogram/generate_triangular_reweighting_temperatures.sh
```

The default output directory is:

```text
data/temperatures/triangular/
```

Both scripts use a temperature spacing of:

```text
0.0005
```

The temperature window itself depends on the lattice size.

## Preparing multihistogram input data

The same input-preparation script is used for both lattice types:

```text
scripts/multihistogram/prepare_eul2d_input.sh
```

It takes three arguments:

```text
1. temperature-files directory
2. simulation-results directory
3. output-data directory
```

### Square lattice

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
  data/temperatures/square \
  results/square_simulations \
  data/eul2d/square
```

### Triangular lattice

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
  data/temperatures/triangular \
  results/triangular_simulations \
  data/eul2d/triangular
```

The script reads simulation outputs from directories of the form:

```text
results/<simulation_type>/T_TTTTT/size_L/
```

For each temperature and lattice size, it combines:

```text
mefile.txt
epfile.txt
enfile.txt
fpfile.txt
fnfile.txt
```

The resulting files are written to:

```text
data/eul2d/<lattice_type>/size_L/conf-T-IDX.dat
```

## Input columns

Each row in `conf-T-IDX.dat` contains ten columns:

1. first value from `mefile.txt`
2. second value from `mefile.txt`
3. value from `epfile.txt`
4. value from `enfile.txt`
5. first value from `fpfile.txt`
6. second value from `fpfile.txt`
7. third value from `fpfile.txt`
8. first value from `fnfile.txt`
9. second value from `fnfile.txt`
10. third value from `fnfile.txt`

The files `bpfile.txt` and `bnfile.txt` produced by the simulations are not used in this input-preparation step.

## Running the multihistogram analysis

The same runner is used for square and triangular lattices:

```text
scripts/multihistogram/run_multihistogram_all.sh
```

It takes:

```text
1. multihistogram executable directory
2. temperature-files directory
3. prepared input-data directory
4. output directory
5. optional refinement factor
```

### Square lattice

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
  src/multihistogram \
  data/temperatures/square \
  data/eul2d/square \
  results/multihistogram/square
```

### Triangular lattice

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
  src/multihistogram \
  data/temperatures/triangular \
  data/eul2d/triangular \
  results/multihistogram/triangular
```

The `Rw` executable receives the prepared input-data directory explicitly, so the same executable can process both lattice types.

## Reweighting grid refinement

The multihistogram runner uses a refinement factor to define the number of output temperature points.

By default:

```text
POINTS = (N_temperatures - 1) * 50 + 1
```

where `N_temperatures` is the number of temperatures in `temperatures_L.txt`.

Each interval between consecutive input temperatures is therefore divided into 50 smaller intervals.

For an input temperature spacing of:

```text
0.0005
```

the corresponding output spacing is approximately:

```text
0.00001
```

A different refinement factor can be passed as the fifth argument.

For example:

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
  src/multihistogram \
  data/temperatures/square \
  data/eul2d/square \
  results/multihistogram/square \
  100
```

## Output files

For each lattice size, the multihistogram runner writes:

```text
results_L.dat
output2_L.dat
```

For example:

```text
results/multihistogram/square/results_128.dat
results/multihistogram/square/output2_128.dat

results/multihistogram/triangular/results_128.dat
results/multihistogram/triangular/output2_128.dat
```

where `L` is the lattice size.

## Generated data

Prepared multihistogram input files and multihistogram outputs are generated locally and are not tracked by Git.

The corresponding large datasets are distributed separately from the source-code repository.
