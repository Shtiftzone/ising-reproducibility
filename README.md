# Ising reproducibility

This repository contains the code and workflow used for two-dimensional Ising model simulations, multihistogram reweighting, and figure generation.

The analysis includes simulations on:

* square lattices
* triangular lattices

The repository is organized as a reproducibility package: source code, portable workflow scripts, documentation, and lightweight data products are kept together, while large generated datasets are distributed separately.

## Repository structure

```text
src/
  simulation/          C simulation code
  multihistogram/      C++ multihistogram reweighting code

scripts/
  simulation/          simulation workflow scripts
  multihistogram/      temperature selection, input preparation, and reweighting scripts

docs/
  simulation_workflow.md
  multihistogram_workflow.md
  software_requirements.md

data/                  generated workflow data directories
results/               generated analysis outputs

requirements.txt       Python dependencies
```

The `data/` and `results/` directories are populated by the workflow and are not part of the source-code repository, except for lightweight data products explicitly included for figure reproduction.

## Workflow overview

The analysis follows the pipeline:

```text
Build simulation executables
        ¡
Run square- or triangular-lattice simulations
        ¡
Generate lattice-specific reweighting temperature files
        ¡
Prepare multihistogram input data
        ¡
Build the multihistogram executable
        ¡
Run multihistogram reweighting
        ¡
Prepare figure data
        ¡
Generate figures
```

Square- and triangular-lattice datasets are kept separate throughout the workflow.

## Quick start

### 1. Build the simulation executables

```bash
cd src/simulation
make
cd ../..
```

### 2. Run the simulations

Square lattice:

```bash
bash scripts/simulation/run_square_2d_simulations.sh \
  src/simulation \
  results/square_simulations \
  2.22 \
  2.32 \
  0.0005 \
  200
```

Triangular lattice:

```bash
bash scripts/simulation/run_triangular_2d_simulations.sh \
  src/simulation \
  results/triangular_simulations \
  3.59 \
  3.72 \
  0.0005 \
  200
```

### 3. Generate the reweighting temperature files

Square lattice:

```bash
bash scripts/multihistogram/generate_square_reweighting_temperatures.sh
```

Triangular lattice:

```bash
bash scripts/multihistogram/generate_triangular_reweighting_temperatures.sh
```

This creates temperature files under:

```text
data/temperatures/square/
data/temperatures/triangular/
```

### 4. Prepare multihistogram input data

Square lattice:

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
  data/temperatures/square \
  results/square_simulations \
  data/eul2d/square
```

Triangular lattice:

```bash
bash scripts/multihistogram/prepare_eul2d_input.sh \
  data/temperatures/triangular \
  results/triangular_simulations \
  data/eul2d/triangular
```

### 5. Build the multihistogram executable

```bash
cd src/multihistogram
make
cd ../..
```

### 6. Run multihistogram reweighting

Square lattice:

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
  src/multihistogram \
  data/temperatures/square \
  data/eul2d/square \
  results/multihistogram/square
```

Triangular lattice:

```bash
bash scripts/multihistogram/run_multihistogram_all.sh \
  src/multihistogram \
  data/temperatures/triangular \
  data/eul2d/triangular \
  results/multihistogram/triangular
```

The default reweighting refinement factor is `50`.

## Documentation

More detailed instructions are available in:

* [`docs/software_requirements.md`](docs/software_requirements.md)
* [`docs/simulation_workflow.md`](docs/simulation_workflow.md)
* [`docs/multihistogram_workflow.md`](docs/multihistogram_workflow.md)

## Python dependencies

Python is used for post-processing and figure generation.

Install the required packages with:

```bash
python -m pip install -r requirements.txt
```

## Data availability

Simulation outputs, prepared multihistogram inputs, and multihistogram results can be large and are generated outside the source-code repository.

Large datasets associated with the analysis are distributed separately from the repository.

Lightweight data products required for reproducing the final figures can be included directly in the repository.

## Citation

Citation information will be provided in `CITATION.cff`.

## License

Licensing information will be provided in the repository `LICENSE` file.
